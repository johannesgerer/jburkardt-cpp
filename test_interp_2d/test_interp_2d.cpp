# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "test_interp_2d.hpp"
# include "r8lib.hpp"

//****************************************************************************80

void f00_f0 ( int fi, int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F00_F0 returns the value of any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FI, the index of the function.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  if ( fi == 1 )
  {
    f01_f0 ( n, x, y, f );
  }
  else if ( fi == 2 )
  {
    f02_f0 ( n, x, y, f );
  }
  else if ( fi == 3 )
  {
    f03_f0 ( n, x, y, f );
  }
  else if ( fi == 4 )
  {
    f04_f0 ( n, x, y, f );
  }
  else if ( fi == 5 )
  {
    f05_f0 ( n, x, y, f );
  }
  else if ( fi == 6 )
  {
    f06_f0 ( n, x, y, f );
  }
  else if ( fi == 7 )
  {
    f07_f0 ( n, x, y, f );
  }
  else if ( fi == 8 )
  {
    f08_f0 ( n, x, y, f );
  }
  else if ( fi == 9 )
  {
    f09_f0 ( n, x, y, f );
  }
  else if ( fi == 10 )
  {
    f10_f0 ( n, x, y, f );
  }
  else if ( fi == 11 )
  {
    f11_f0 ( n, x, y, f );
  }
  else if ( fi == 12 )
  {
    f12_f0 ( n, x, y, f );
  }
  else if ( fi == 13 )
  {
    f13_f0 ( n, x, y, f );
  }
  else
  {
    cerr << "\n";
    cerr << "F00_F0 - Fatal error!\n";
    cerr << "  Illegal function index FI = " << fi << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void f00_f1 ( int fi, int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F00_F1 returns first derivatives of any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FI, the index of the function.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the first derivative values.
//
{
  if ( fi == 1 )
  {
    f01_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 2 )
  {
    f02_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 3 )
  {
    f03_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 4 )
  {
    f04_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 5 )
  {
    f05_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 6 )
  {
    f06_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 7 )
  {
    f07_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 8 )
  {
    f08_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 9 )
  {
    f09_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 10 )
  {
    f10_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 11 )
  {
    f11_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 12 )
  {
    f12_f1 ( n, x, y, fx, fy );
  }
  else if ( fi == 13 )
  {
    f13_f1 ( n, x, y, fx, fy );
  }
  else
  {
    cerr << "\n";
    cerr << "F00_F1 - Fatal error!\n";
    cerr << "  Illegal function index FI = " << fi << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void f00_f2 ( int fi, int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F00_F2 returns second derivatives of any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FI, the index of the function.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], the second derivatives.
//
{
  if ( fi == 1 )
  {
    f01_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 2 )
  {
    f02_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 3 )
  {
    f03_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 4 )
  {
    f04_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 5 )
  {
    f05_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 6 )
  {
    f06_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 7 )
  {
    f07_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 8 )
  {
    f08_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 9 )
  {
    f09_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 10 )
  {
    f10_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 11 )
  {
    f11_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 12 )
  {
    f12_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else if ( fi == 13 )
  {
    f13_f2 ( n, x, y, fxx, fxy, fyy );
  }
  else
  {
    cerr << "\n";
    cerr << "F00_F2 - Fatal error!\n";
    cerr << "  Illegal function index FI = " << fi << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

int f00_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    F00_NUM returns the number of test functions available.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//   Output, int FUN_NUM, the number of test functions.
//
{
  int fun_num;

  fun_num = 13;

  return fun_num;
}
//****************************************************************************80

string f00_title ( int fi )

//****************************************************************************80
//
//  Purpose:
//
//    F00_TITLE returns the title for any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FI, the index of the function.
//
//    Output, string FT, the function title.
//
{
  string ft;

  if ( fi == 1 )
  {
    ft = f01_title ( );
  }
  else if ( fi == 2 )
  {
    ft = f02_title ( );
  }
  else if ( fi == 3 )
  {
    ft = f03_title ( );
  }
  else if ( fi == 4 )
  {
    ft = f04_title ( );
  }
  else if ( fi == 5 )
  {
    ft = f05_title ( );
  }
  else if ( fi == 6 )
  {
    ft = f06_title ( );
  }
  else if ( fi == 7 )
  {
    ft = f07_title ( );
  }
  else if ( fi == 8 )
  {
    ft = f08_title ( );
  }
  else if ( fi == 9 )
  {
    ft = f09_title ( );
  }
  else if ( fi == 10 )
  {
    ft = f10_title ( );
  }
  else if ( fi == 11 )
  {
    ft = f11_title ( );
  }
  else if ( fi == 12 )
  {
    ft = f12_title ( );
  }
  else if ( fi == 13 )
  {
    ft = f13_title ( );
  }
  else
  {
    cerr << "\n";
    cerr << "F00_TITLE - Fatal error!\n";
    cerr << "  Illegal function index FI = " << fi << "\n";
    exit ( 1 );
  }

  return ft;
}
//****************************************************************************80

void f01_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F01_F0 returns the value of function 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = 
      0.75 * exp ( - ( pow ( 9.0 * x[i] - 2.0, 2 )            
                     + pow ( 9.0 * y[i] - 2.0, 2 ) ) / 4.0 )  
    + 0.75 * exp ( - ( pow ( 9.0 * x[i] + 1.0, 2 ) ) / 49.0   
                     - ( 9.0 * y[i] + 1.0 ) / 10.0 )      
    + 0.5  * exp ( - ( pow ( 9.0 * x[i] - 7.0, 2 )            
                     + pow ( 9.0 * y[i] - 3.0, 2 ) ) / 4.0 )  
    - 0.2  * exp (   - pow ( 9.0 * x[i] - 4.0, 2 )            
                     - pow ( 9.0 * y[i] - 7.0, 2 ) );
  }
  return;
}
//****************************************************************************80

void f01_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F01_F1 returns first derivatives of function 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4;  

  for ( i = 0; i < n; i++ )
  {
    t1 = exp ( - ( pow ( 9.0 * x[i] - 2.0, 2 )
                 + pow ( 9.0 * y[i] - 2.0, 2 ) ) / 4.0 );
    t2 = exp ( - ( pow ( 9.0 * x[i] + 1.0, 2 ) ) / 49.0
                     - ( 9.0 * y[i] + 1.0 ) / 10.0 );
    t3 = exp ( - ( pow ( 9.0 * x[i] - 7.0, 2 )
                 + pow ( 9.0 * y[i] - 3.0, 2 ) ) / 4.0 );
    t4 = exp ( -   pow ( 9.0 * x[i] - 4.0, 2 )
                 - pow ( 9.0 * y[i] - 7.0, 2 ) );

    fx[i] =
     - 3.375           * ( 9.0 * x[i] - 2.0 ) * t1
     - ( 27.0 / 98.0 ) * ( 9.0 * x[i] + 1.0 ) * t2
     - 2.25            * ( 9.0 * x[i] - 7.0 ) * t3
     + 3.6             * ( 9.0 * x[i] - 4.0 ) * t4;

    fy[i] =
     - 3.375 * ( 9.0 * y[i] - 2.0 ) * t1
     - 0.675                        * t2
     - 2.25  * ( 9.0 * y[i] - 3.0 ) * t3
     + 3.6   * ( 9.0 * y[i] - 7.0 ) * t4;

  }
  return;
}
//****************************************************************************80

void f01_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F01_F2 returns second derivatives of function 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4;

  for ( i = 0; i < n; i++ )
  {
    t1 = exp ( - ( pow ( 9.0 * x[i] - 2.0, 2 )
                 + pow ( 9.0 * y[i] - 2.0, 2 ) ) / 4.0 );
    t2 = exp ( - ( pow ( 9.0 * x[i] + 1.0, 2 ) ) / 49.0
                     - ( 9.0 * y[i] + 1.0 ) / 10.0 );
    t3 = exp ( - ( pow ( 9.0 * x[i] - 7.0, 2 )
                 + pow ( 9.0 * y[i] - 3.0, 2 ) ) / 4.0 );
    t4 = exp ( -   pow ( 9.0 * x[i] - 4.0, 2 )
                 - pow ( 9.0 * y[i] - 7.0, 2 ) );

    fxx[i] =
       15.1875 * ( pow ( 9.0 * x[i] - 2.0, 2 ) - 2.0 )  * t1
     + 60.75   * ( pow ( 9.0 * x[i] + 1.0, 2 ) - 24.5 ) * t2
     + 10.125  * ( pow ( 9.0 * x[i] - 7.0, 2 ) - 2.0 )  * t3
     - 64.8    * ( pow ( 9.0 * x[i] - 4.0, 2 ) - 0.5 )  * t4;

    fxy[i] =
       15.1875 * ( 9.0 * x[i] - 2.0 ) * ( 9.0 * y[i] - 2.0 ) * t1
     + ( 243.0 / 980.0 ) * ( 9.0 * x[i] + 1.0 ) * t2
     + 10.125 * ( 9.0 * x[i] - 7.0 ) * ( 9.0 * y[i] - 3.0 ) * t3
     - 64.8 * ( 9.0 * x[i] - 4.0 ) * ( 9.0 * y[i] - 7.0 ) * t4;

    fyy[i] =
       15.1875 * ( pow ( 9.0 * y[i] - 2.0, 2 ) - 2.0 ) * t1
     + 0.6075  *                                     t2
     + 10.125  * ( pow ( 9.0 * y[i] - 3.0, 2 ) - 2.0 ) * t3
     - 64.8    * ( pow ( 9.0 * y[i] - 7.0, 2 ) - 0.5 ) * t4;

  }
  return;
}
//****************************************************************************80

string f01_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F01_TITLE returns the title for function 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Exponential";

  return ft;
}
//****************************************************************************80

void f02_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F02_F0 returns the value of function 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = ( tanh ( 9.0 * ( y[i] - x[i] ) ) + 1.0 ) / 9.0;
  }

  return;
}
//****************************************************************************80

void f02_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F02_F1 returns first derivatives of function 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i;
  double t1;

  for ( i = 0; i < n; i++ )
  {
    t1 = 18.0 * ( y[i] - x[i] );
    fx[i] = - 4.0 / ( exp ( t1 ) + 2.0 + exp ( - t1 ) );
    fy[i] = - fx[i];
  }
  return;
}
//****************************************************************************80

void f02_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F02_F2 returns second derivatives of function 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;
  double t1;

  for ( i = 0; i < n; i++ )
  {
    t1 = 18.0 * ( y[i] - x[i] );

    fxx[i] = 18.0 * tanh ( 0.5 * t1 ) 
      * ( tanh ( 9.0 * ( y[i] - x[i] ) ) + 1.0 ) / 9.0;
    fxy[i] = - fxx[i];
    fyy[i] = fxx[i];
  }
  return;
}
//****************************************************************************80

string f02_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F02_TITLE returns the title for function 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Cliff";

  return ft;
}
//****************************************************************************80

void f03_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F03_F0 returns the value of function 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = ( 1.25 + cos ( 5.4 * y[i] ) )
     / ( 6.0 + 6.0 * pow ( 3.0 * x[i] - 1.0, 2 ) );
  }

  return;
}
//****************************************************************************80

void f03_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F03_F1 returns first derivatives of function 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i;
  double t1;
  double t2;  

  for ( i = 0; i < n; i++ )
  {
    t1 = 5.4 * y[i];
    t2 = 1.0 + pow ( 3.0 * x[i] - 1.0, 2 );
    fx[i] = - ( 3.0 * x[i] - 1.0 ) * ( 1.25 + cos ( t1 ) ) / ( t2 * t2 );
    fy[i] = - 0.9 * sin ( t1 ) / t2;
  }
  return;
}
//****************************************************************************80

void f03_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F03_F2 returns second derivatives of function 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;
  double t1;
  double t2; 

  for ( i = 0; i < n; i++ ) 
  {
    t1 = 5.4 * y[i];
    t2 = 1.0 + pow ( 3.0 * x[i] - 1.0, 2 );

    fxx[i] = 3.0 * ( 1.25 + cos ( t1 ) ) * ( 3.0 * t2 - 4.0 ) / pow ( t2, 3 );
    fxy[i] = 5.4 * ( 3.0 * x[i] - 1.0 ) * sin ( t1 ) / ( t2 * t2 );
    fyy[i] = - 4.86 * cos ( t1 ) / t2;
  }
  return;
}
//****************************************************************************80

string f03_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F03_TITLE returns the title for function 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Saddle";

  return ft;
}
//****************************************************************************80

void f04_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F04_F0 returns the value of function 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = exp ( - 5.0625 * ( pow ( x[i] - 0.5, 2 )
                            + pow ( y[i] - 0.5, 2 ) ) ) / 3.0;
  }
  return;
}
//****************************************************************************80

void f04_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F04_F1 returns first derivatives of function 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i;
  double t1;
  double t2;
  double t3;

  for ( i = 0; i < n; i++ )
  {
    t1 = x[i] - 0.5;
    t2 = y[i] - 0.5;
    t3 = - 3.375 * exp ( - 5.0625 * ( t1 * t1 + t2 * t2 ) );
    fx[i] = t1 * t3;
    fy[i] = t2 * t3;
  }

  return;
}
//****************************************************************************80

void f04_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F04_F2 returns second derivatives of function 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;
  double t1;
  double t2;
  double t3;

  for ( i = 0; i < n; i++ )
  {
    t1 = x[i] - 0.5;
    t2 = y[i] - 0.5;
    t3 = - 3.375 * exp ( - 5.0625 * ( t1 * t1 + t2 * t2 ) );

    fxx[i] = ( 1.0 - 10.125 * t1 * t1 ) * t3;
    fxy[i] = - 10.125 * t1 * t2 * t3;
    fyy[i] = ( 1.0 - 10.125 * t2 * t2 ) * t3;
  }
  return;
}
//****************************************************************************80

string f04_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F04_TITLE returns the title for function 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Gentle";

  return ft;
}
//****************************************************************************80

void f05_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F05_F0 returns the value of function 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i; 

  for ( i = 0; i < n; i++ )
  {
    f[i] = exp ( - 20.25 * ( pow ( x[i] - 0.5, 2 )
                           + pow ( y[i] - 0.5, 2 ) ) ) / 3.0;
  }
  return;
}
//****************************************************************************80

void f05_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F05_F1 returns first derivatives of function 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i;
  double t1;
  double t2;
  double t3;

  for ( i = 0; i < n; i++ )
  {
    t1 = x[i] - 0.5;
    t2 = y[i] - 0.5;
    t3 = - 13.5 * exp ( - 20.25 * ( t1 * t1 + t2 * t2 ) );
    fx[i] = t1 * t3;
    fy[i] = t2 * t3;
  }

  return;
}
//****************************************************************************80

void f05_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F05_F2 returns second derivatives of function 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;
  double t1;
  double t2;
  double t3;  

  for ( i = 0; i < n; i++ )
  {
    t1 = x[i] - 0.5;
    t2 = y[i] - 0.5;
    t3 = - 13.5 * exp ( - 20.25 * ( t1 * t1 + t2 * t2 ) );

    fxx[i] = ( 1.0 - 40.5 * t1 * t1 ) * t3;
    fxy[i] = - 40.5 * t1 * t2 * t3;
    fyy[i] = ( 1.0 - 40.5 * t2 * t2 ) * t3;
  }
  return;
}
//****************************************************************************80

string f05_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F05_TITLE returns the title for function 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Steep";

  return ft;
}
//****************************************************************************80

void f06_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F06_F0 returns the value of function 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;
  double t4;  

  for ( i = 0; i < n; i++ )
  {
    t4 = 64.0 - 81.0 * ( pow ( x[i] - 0.5, 2 ) + pow ( y[i] - 0.5, 2 ) );

    if ( 0.0 <= t4 )
    {
      f[i] = sqrt ( t4 ) / 9.0 - 0.5;
    }
    else
    {
      f[i] = 0.0;
    }
  }
  return;
}
//****************************************************************************80

void f06_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F06_F1 returns first derivatives of function 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4;

  for ( i = 0; i < n; i++ )
  {
    t4 = 64.0 - 81.0 * ( pow ( x[i] - 0.5, 2 ) + pow ( y[i] - 0.5, 2 ) );

    if ( 0.0 < t4 )
    {
      t1 = x[i] - 0.5;
      t2 = y[i] - 0.5;
      t3 = - 9.0 / sqrt ( t4 );
      fx[i] = t1 * t3;
      fy[i] = t2 * t3;
    }
    else
    {
      fx[i] = 0.0;
      fy[i] = 0.0;
    }
  }
  return;
}
//****************************************************************************80

void f06_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F06_F2 returns second derivatives of function 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4;

  for ( i = 0; i < n; i++ )
  {
    t4 = 64.0 - 81.0 * ( pow ( x[i] - 0.5, 2 ) + pow ( y[i] - 0.5, 2 ) );

    if ( 0.0 < t4 )
    {
      t1 = x[i] - 0.5;
      t2 = y[i] - 0.5;
      t3 = - 9.0 / sqrt ( t4 );
      fxx[i] = ( 1.0 + t1 * t3 * t1 * t3 ) * t3;
      fxy[i] = t1 * t3 * t2 * t3;
      fyy[i] = ( 1.0 + t2 * t3 * t2 * t3 ) * t3;
    }
    else
    {
      fxx[i] = 0.0;
      fxy[i] = 0.0;
      fyy[i] = 0.0;
    }
  }
  return;
}
//****************************************************************************80

string f06_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F06_TITLE returns the title for function 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Sphere";

  return ft;
}
//****************************************************************************80

void f07_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F07_F0 returns the value of function 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;  

  for ( i = 0; i < n; i++ )
  {
    f[i] = 2.0 * cos ( 10.0 * x[i] ) * sin ( 10.0 * y[i] )
     + sin ( 10.0 * x[i] * y[i] );
  }
  return;
}
//****************************************************************************80

void f07_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F07_F1 returns first derivatives of function 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i;
  double t1;
  double t2;
  double t3;

  for ( i = 0; i < n; i++ )
  {
    t1 = 10.0 * x[i];
    t2 = 10.0 * y[i];
    t3 = 10.0 * cos ( 10.0 * x[i] * y[i] );
    fx[i] = - 20.0 * sin ( t1 ) * sin ( t2 ) + t3 * y[i];
    fy[i] =   20.0 * cos ( t1 ) * cos ( t2 ) + t3 * x[i];
  }
  return;
}
//****************************************************************************80

void f07_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F07_F2 returns second derivatives of function 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4;  

  for ( i = 0; i < n; i++ )
  {
    t1 = 10.0 * x[i];
    t2 = 10.0 * y[i];
    t3 = 10.0 * cos ( 10.0 * x[i] * y[i] );
    t4 = 100.0 * sin ( 10.0 * x[i] * y[i] );

    fxx[i] = - 200.0 * cos ( t1 ) * sin ( t2 )      - t4 * y[i] * y[i];
    fxy[i] = - 200.0 * sin ( t1 ) * cos ( t2 ) + t3 - t4 * x[i] * y[i];
    fyy[i] = - 200.0 * cos ( t1 ) * sin ( t2 )      - t4 * x[i] * x[i];
  }
  return;
}
//****************************************************************************80

string f07_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F07_TITLE returns the title for function 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Trig";

  return ft;
}
//****************************************************************************80

void f08_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F08_F0 returns the value of function 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4; 

  for ( i = 0; i < n; i++ )
  {
    t1 = 5.0 - 10.0 * x[i];
    t2 = 5.0 - 10.0 * y[i];
    t3 = exp ( - 0.5 * t1 * t1 );
    t4 = exp ( - 0.5 * t2 * t2 );
    f[i] = t3 + 0.75 * t4 * ( 1.0 + t3 );
  }
  return;
}
//****************************************************************************80

void f08_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F08_F1 returns first derivatives of function 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4;  

  for ( i = 0; i < n; i++ )
  {
    t1 = 5.0 - 10.0 * x[i];
    t2 = 5.0 - 10.0 * y[i];
    t3 = exp ( - 0.5 * t1 * t1 );
    t4 = exp ( - 0.5 * t2 * t2 );

    fx[i] = t1 * t3 * ( 10.0 + 7.5 * t4 );
    fy[i] = t2 * t4 * (  7.5 + 7.5 * t3 );
  }
  return;
}
//****************************************************************************80

void f08_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F08_F2 returns second derivatives of function 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4;

  for ( i = 0; i < n; i++ )
  {
    t1 = 5.0 - 10.0 * x[i];
    t2 = 5.0 - 10.0 * y[i];
    t3 = exp ( - 0.5 * t1 * t1 );
    t4 = exp ( - 0.5 * t2 * t2 );

    fxx[i] = t3 * ( t1 * t1 - 1.0 ) * ( 100.0 + 75.0 * t4 );
    fxy[i] = 75.0 * t1 * t2 * t3 * t4;
    fyy[i] = t4 * ( t2 * t2 - 1.0 ) * ( 75.0 + 75.0 * t3 );
  }
  return;
}
//****************************************************************************80

string f08_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F08_TITLE returns the title for function 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Synergistic Gaussian";

  return ft;
}
//****************************************************************************80

void f09_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F09_F0 returns the value of function 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4; 

  for ( i = 0; i < n; i++ )
  {
    t1 = exp ( ( 10.0 - 20.0 * x[i] ) / 3.0 );
    t2 = exp ( ( 10.0 - 20.0 * y[i] ) / 3.0 );
    t3 = 1.0 / ( 1.0 + t1 );
    t4 = 1.0 / ( 1.0 + t2 );
    f[i] = pow ( pow ( 20.0 / 3.0, 3 ) * t1 * t2, 2 )
     * pow ( t3 * t4, 5 )
     * ( t1 - 2.0 * t3 ) * ( t2 - 2.0 * t4 );
  }
  return;
}
//****************************************************************************80

void f09_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F09_F1 returns first derivatives of function 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4;

  for ( i = 0; i < n; i++ )
  {
    t1 = exp ( ( 10.0 - 20.0 * x[i] ) / 3.0 );
    t2 = exp ( ( 10.0 - 20.0 * y[i] ) / 3.0 );
    t3 = 1.0 / ( 1.0 + t1 );
    t4 = 1.0 / ( 1.0 + t2 );

    fx[i] = pow ( ( 20.0 / 3.0 ) * t1, 2 ) * pow ( ( 20.0 / 3.0 ) * t3, 5 )
     * ( 2.0 * t1 - 3.0 * t3 - 5.0 + 12.0 * t3 * t3 )
     * t2 * t2 * pow ( t4, 5 )
     * ( t2 - 2.0 * t4 );

    fy[i] = pow ( ( 20.0 / 3.0 ) * t1, 2 ) * pow ( ( 20.0 / 3.0 ) * t3, 5 )
     * ( 2.0 * t2 - 3.0 * t4 - 5.0 + 12.0 * t4 * t4 )
     * t2 * t2 * pow ( t4, 5 )
     * ( t1 - 2.0 * t3 );
  }
  return;
}
//****************************************************************************80

void f09_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F09_F2 returns second derivatives of function 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;  

  for ( i = 0; i < n; i++ )
  {
    t1 = exp ( ( 10.0 - 20.0 * x[i] ) / 3.0 );
    t2 = exp ( ( 10.0 - 20.0 * y[i] ) / 3.0 );
    t3 = 1.0 / ( 1.0 + t1 );
    t4 = 1.0 / ( 1.0 + t2 );
    t5 = 20.0 / 3.0;
    t6 = pow ( t5 * t1 * t2, 2 ) * pow ( t5 * t3 * t4, 5 );

    fxx[i] = t5 * t6 * ( t2 - 2.0 * t4 )
     * ( ( ( - 84.0 * t3 + 78.0 ) * t3 + 23.0 ) * t3
     + 4.0 * t1 - 25.0 );

    fxy[i] = t5 * t6
     * ( ( 12.0 * t4 - 3.0 ) * t4 + 2.0 * t2 - 5.0 )
     * ( ( 12.0 * t3 - 3.0 ) * t3 + 2.0 * t1 - 5.0 );

    fyy[i] = t5 * t6 * ( t1 - 2.0 * t3 )
     * ( ( ( - 84.0 * t4 + 78.0 ) * t4 + 23.0 ) * t4
     + 4.0 * t2 - 25.0 );
  }
  return;
}
//****************************************************************************80

string f09_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F09_TITLE returns the title for function 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Cloverleaf Asymmetric Peak/Valley";

  return ft;
}
//****************************************************************************80

void f10_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F10_F0 returns the value of function f10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;
  double t1;
  double t2;
  double t3;  

  for ( i = 0; i < n; i++ )
  {
    t1 = sqrt ( pow ( 80.0 * x[i] - 40.0, 2 ) 
              + pow ( 90.0 * y[i] - 45.0, 2 ) );
    t2 = exp ( - 0.04 * t1 );
    t3 = cos ( 0.15 * t1 );

    f[i] = t2 * t3;
  }
  return;
}
//****************************************************************************80

void f10_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F10_F1 returns first derivatives of function f10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4; 

  for ( i = 0; i < n; i++ )
  {
    t1 = sqrt ( pow ( 80.0 * x[i] - 40.0, 2 )
              + pow ( 90.0 * y[i] - 45.0, 2 ) );

    if ( t1 == 0.0 )
    {
      fx[i] = 0.0;
      fy[i] = 0.0;
    }
    else
    {
      t2 = exp ( - 0.04 * t1 );
      t3 = cos ( 0.15 * t1 );
      t4 = sin ( 0.15 * t1 );
      fx[i] = - t2 * ( 12.0 * t4 + 3.2 * t3 ) 
       * ( 80.0 * x[i] - 40.0 ) / t1;
      fy[i] = - t2 * ( 13.5 * t4 + 3.6 * t3 ) 
       * ( 90.0 * y[i] - 45.0 ) / t1;
    }
  }
  return;
}
//****************************************************************************80

void f10_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F10_F2 returns second derivatives of function f10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;

  for ( i = 0; i < n; i++ )
  {
    t1 = sqrt ( pow ( 80.0 * x[i] - 40.0, 2 )
              + pow ( 90.0 * y[i] - 45.0, 2 ) );

    if ( t1 == 0.0 )
    {
      fxx[i] = 0.0;
      fxy[i] = 0.0;
      fyy[i] = 0.0;
    }
    else
    {
      t2 = exp ( - 0.04 * t1 );
      t3 = cos ( 0.15 * t1 );
      t4 = sin ( 0.15 * t1 );
      t5 = t2 / pow ( t1, 3 );

      fxx[i] = t5 * ( t1 * ( 76.8 * t4 - 133.76 * t3 )
       * pow ( 80.0 * x[i] - 40.0, 2 )
       - ( 960.0 * t4 + 256.0 * t3 ) 
       * pow ( 90.0 * y[i] - 45.0, 2 ) );

      fxy[i] = t5 * ( t1 * ( 86.4 * t4 - 150.48 * t3 ) 
       + 1080.0 * t4 + 288.0 * t3 ) * ( 80.0 * x[i] - 40.0 ) 
       * ( 90.0 * y[i] - 45.0 );

      fyy[i] = t5 * ( t1 * ( 97.2 * t4 - 169.29 * t3 )
       * pow ( 90.0 * y[i] - 45.0, 2 ) - ( 1215.0 * t4 + 324.0 * t3 ) 
       * pow ( 80.0 * x[i] - 40.0, 2 ) );
    }
  }
  return;
}
//****************************************************************************80

string f10_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F10_TITLE returns the title for function f10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Cosine Peak";

  return ft;
}
//****************************************************************************80

void f11_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F11_F0 returns the value of function f11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = x[i] * ( y[i] + 1.0 );
  }
  return;
}
//****************************************************************************80

void f11_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F11_F1 returns first derivatives of function f11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i; 

  for ( i = 0; i < n; i++ )
  {
    fx[i] = y[i] + 1.0;
    fy[i] = x[i];
  }
  return;
}
//****************************************************************************80

void f11_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F11_F2 returns second derivatives of function f11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fxx[i] = 0.0;
    fxy[i] = 1.0;
    fyy[i] = 0.0;
  }
  return;
}
//****************************************************************************80

string f11_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F11_TITLE returns the title for function f11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Bilinear function";

  return ft;
}
//****************************************************************************80

void f12_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F12_F0 returns the value of function f12.
//
//  Discussion:
//
//    This is an example from Vicente Romero.
//
//      R = sqrt ( X^2 + Y^2 )
//      T = atan ( Y / X )
//      F(X,Y) = ( 0.8 * R + 0.35 * sin ( 2.4 * pi * R / sqrt ( 2 )  ) )
//               * 1.5 * sin ( 1.3 * T )
//
//    The mean and standard deviation of the function over the interval
//    are approximately:
//
//      mu    = 0.581608
//      sigma = 0.343208
//
//    Since the interpolation interval is the unit square, this means the
//    integral of the function over the interval can also be estimated as
//
//      I = 0.581608
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Vicente Romero, John Burkardt, Max Gunzburger, Janet Peterson,
//    Initial Evaluation of Centroidal Voronoi Tessellation for
//    Statistical Sampling and Function Integration,
//    Fourth International Symposium on Uncertainty Modeling and Analysis,
//    (ISUMA) 2003.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;
  static double pi = 3.141592653589793;
  double r;
  double t;

  for ( i = 0; i < n; i++ )
  {
    r = sqrt ( x[i] * x[i] + y[i] * y[i] );
    t = atan2 ( y[i], x[i] );

    f[i] = 1.5 * ( 0.8 * r 
    + 0.35 * sin ( 2.4 * pi * r / sqrt ( 2.0 ) ) ) 
    * sin ( 1.3 * t );
  }
  return;
}
//****************************************************************************80

void f12_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F12_F1 returns first derivatives of function f12.
//
//  Discussion:
//
//    Currently, the derivative information is of no interest to me.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i; 

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 0.0;
    fy[i] = 0.0;
  }
  return;
}
//****************************************************************************80

void f12_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F12_F2 returns second derivatives of function f12.
//
//  Discussion:
//
//    Currently, the derivative information is of no interest to me.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fxx[i] = 0.0;
    fxy[i] = 0.0;
    fyy[i] = 0.0;
  }
  return;
}
//****************************************************************************80

string f12_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F12_TITLE returns the title for function f12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Vicente Romero function";

  return ft;
}
//****************************************************************************80

void f13_f0 ( int n, double x[], double y[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    F13_F0 returns the value of function f13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double F[N], the function values.
//
{
  int i;
  static double pi = 3.141592653589793;
  double r;
  double t;

  for ( i = 0; i < n; i++ )
  {
    f[i] = 1.0 / ( pow ( 10.0 * x[i] - 5.0, 2 ) 
                 + pow ( 10.0 * y[i] - 5.0, 2 ) + 1.0 );
  }
  return;
}
//****************************************************************************80

void f13_f1 ( int n, double x[], double y[], double fx[], double fy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F13_F1 returns first derivatives of function f13.
//
//  Discussion:
//
//    Currently, the derivative information is of no interest to me.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FX[N], FY[N], the derivative values.
//
{
  int i; 

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 0.0;
    fy[i] = 0.0;
  }
  return;
}
//****************************************************************************80

void f13_f2 ( int n, double x[], double y[], double fxx[], double fxy[], 
  double fyy[] )

//****************************************************************************80
//
//  Purpose:
//
//    F13_F2 returns second derivatives of function f13.
//
//  Discussion:
//
//    Currently, the derivative information is of no interest to me.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the evalution points.
//
//    Output, double FXX[N], FXY[N], FYY[N], second derivatives.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    fxx[i] = 0.0;
    fxy[i] = 0.0;
    fyy[i] = 0.0;
  }
  return;
}
//****************************************************************************80

string f13_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    F13_TITLE returns the title for function f13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string FT, the function title.
//
{
  string ft;

  ft = "Rescaled Runge function";

  return ft;
}
//****************************************************************************80

int g00_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    G00_NUM returns the number of grids available.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//   Output, int GRID_NUM, the number of grids.
//
{
  int grid_num;

  grid_num = 5;

  return grid_num;
}
//****************************************************************************80

int g00_size ( int gi )

//****************************************************************************80
//
//  Purpose:
//
//    G00_SIZE returns the size for any grid.
//
//  Discussion:
//
//    The "grid size" is simply the number of points in the grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int GI, the index of the grid.
//
//    Output, int GN, the grid size.
//
{
  int gn;

  if ( gi == 1 )
  {
    gn = g01_size ( );
  }
  else if ( gi == 2 )
  {
    gn = g02_size ( );
  }
  else if ( gi == 3 )
  {
    gn = g03_size ( );
  }
  else if ( gi == 4 )
  {
    gn = g04_size ( );
  }
  else if ( gi == 5 )
  {
    gn = g05_size ( );
  }
  else
  {
    cout << "\n";
    cout << "G00_SIZE - Fatal error!\n";
    cout << "  Illegal grid index GI = " << gi << "\n";
    exit ( 1 );
  }
  return gn;
}
//****************************************************************************80

string g00_title ( int gi )

//****************************************************************************80
//
//  Purpose:
//
//    G00_TITLE returns the title for any grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int GI, the index of the grid.
//
//    Output, string GT, the grid title.
//
{
  string gt;

  if ( gi == 1 )
  {
    gt = g01_title ( );
  }
  else if ( gi == 2 )
  {
    gt = g02_title ( );
  }
  else if ( gi == 3 )
  {
    gt = g03_title ( );
  }
  else if ( gi == 4 )
  {
    gt = g04_title ( );
  }
  else if ( gi == 5 )
  {
    gt = g05_title ( );
  }
  else
  {
    cout << "\n";
    cout << "G00_TITLE - Fatal error!\n";
    cout << "  Illegal grid index GI = " << gi << "\n";
    exit ( 1 );
  }
  return gt;
}
//****************************************************************************80

void g00_xy ( int gi, int gn, double gx[], double gy[] )

//****************************************************************************80
//
//  Purpose:
//
//    G00_XY returns the grid points for any grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int GI, the index of the grid.
//
//    Input, int GN, the grid size.
//
//    Output, double GX[GN], GY[GN], the grid coordinates.
//
{
  if ( gi == 1 )
  {
    g01_xy ( gn, gx, gy );
  }
  else if ( gi == 2 )
  {
    g02_xy ( gn, gx, gy );
  }
  else if ( gi == 3 )
  {
    g03_xy ( gn, gx, gy );
  }
  else if ( gi == 4 )
  {
    g04_xy ( gn, gx, gy );
  }
  else if ( gi == 5 )
  {
    g05_xy ( gn, gx, gy );
  }
  else
  {
    cout << "\n";
    cout << "G00_XY - Fatal error!\n";
    cout << "  Illegal grid index GI = " << gi << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80
 
int g01_size ( )

//****************************************************************************80
//
//  Purpose:
//
//    G01_SIZE returns the size for grid 1.
//
//  Discussion:
//
//    The "grid size" is simply the number of points in the grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GN, the grid size.
//
{
  int gn;

  gn = 100;

  return gn;
}
//****************************************************************************80

string g01_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    G01_TITLE returns the title for grid 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string GT, the grid title.
//
{
  string gt;

  gt = "Franke's 100 node set";

  return gt;
}
//****************************************************************************80

void g01_xy ( int gn, double gx[], double gy[] )

//****************************************************************************80
//
//  Purpose:
//
//    G01_XY returns the grid points for grid 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int GN, the grid size.
//
//    Output, double GX[GN], GY[GN], the grid coordinates.
//
{
  double x[100] = {
    0.0227035,  0.0539888,  0.0217008,  0.0175129,  0.0019029,
   -0.0509685,  0.0395408, -0.0487061,  0.0315828, -0.0418785,
    0.1324189,  0.1090271,  0.1254439,  0.0934540,  0.0767578,
    0.1451874,  0.0626494,  0.1452734,  0.0958668,  0.0695559,
    0.2645602,  0.2391645,  0.2088990,  0.2767329,  0.1714726,
    0.2266781,  0.1909212,  0.1867647,  0.2304634,  0.2426219,
    0.3663168,  0.3857662,  0.3832392,  0.3179087,  0.3466321,
    0.3776591,  0.3873159,  0.3812917,  0.3795364,  0.2803515,
    0.4149771,  0.4277679,  0.4200010,  0.4663631,  0.4855658,
    0.4092026,  0.4792578,  0.4812279,  0.3977761,  0.4027321,
    0.5848691,  0.5730076,  0.6063893,  0.5013894,  0.5741311,
    0.6106955,  0.5990105,  0.5380621,  0.6096967,  0.5026188,
    0.6616928,  0.6427836,  0.6396475,  0.6703963,  0.7001181,
    0.6333590,  0.6908947,  0.6895638,  0.6718889,  0.6837675,
    0.7736939,  0.7635332,  0.7410424,  0.8258981,  0.7306034,
    0.8086609,  0.8214531,  0.7290640,  0.8076643,  0.8170951,
    0.8424572,  0.8684053,  0.8366923,  0.9418461,  0.8478122,
    0.8599583,  0.9175700,  0.8596328,  0.9279871,  0.8512805,
    1.0449820,  0.9670631,  0.9857884,  0.9676313,  1.0129299,
    0.9657040,  1.0019855,  1.0359297,  1.0414677,  0.9471506 };
  double y[100] = {
   -0.0310206,   0.1586742,   0.2576924,   0.3414014,   0.4943596,
    0.5782854,   0.6993418,   0.7470194,   0.9107649,   0.9962890,
    0.0501330,   0.0918555,   0.2592973,   0.3381592,   0.4171125,
    0.5615563,   0.6552235,   0.7524066,   0.9146523,   0.9632421,
    0.0292939,   0.0602303,   0.2668783,   0.3696044,   0.4801738,
    0.5940595,   0.6878797,   0.8185576,   0.9046507,   0.9805412,
    0.0396955,   0.0684484,   0.2389548,   0.3124129,   0.4902989,
    0.5199303,   0.6445227,   0.8203789,   0.8938079,   0.9711719,
   -0.0284618,   0.1560965,   0.2262471,   0.3175094,   0.3891417,
    0.5084949,   0.6324247,   0.7511007,   0.8489712,   0.9978728,
   -0.0271948,   0.1272430,   0.2709269,   0.3477728,   0.4259422,
    0.6084711,   0.6733781,   0.7235242,   0.9242411,   1.0308762,
    0.0255959,   0.0707835,   0.2008336,   0.3259843,   0.4890704,
    0.5096324,   0.6697880,   0.7759569,   0.9366096,   1.0064516,
    0.0285374,   0.1021403,   0.1936581,   0.3235775,   0.4714228,
    0.6091595,   0.6685053,   0.8022808,   0.8476790,   1.0512371,
    0.0380499,   0.0902048,   0.2083092,   0.3318491,   0.4335632,
    0.5910139,   0.6307383,   0.8144841,   0.9042310,   0.9696030,
   -0.0120900,   0.1334114,   0.2695844,   0.3795281,   0.4396054,
    0.5044425,   0.6941519,   0.7459923,   0.8682081,   0.9801409 };

  r8vec_copy ( gn, x, gx );
  r8vec_copy ( gn, y, gy );

  return;
}
//****************************************************************************80

int g02_size ( )

//****************************************************************************80
//
//  Purpose:
//
//    G02_SIZE returns the size for grid 2.
//
//  Discussion:
//
//    The "grid size" is simply the number of points in the grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GN, the grid size.
//
{
  int gn;

  gn = 33;

  return gn;
}
//****************************************************************************80

string g02_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    G02_TITLE returns the title for grid 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string GT, the grid title.
//
{
  string gt;

  gt = "Franke's 33 node set";

  return gt;
}
//****************************************************************************80

void g02_xy ( int gn, double gx[], double gy[] )

//****************************************************************************80
//
//  Purpose:
//
//    G02_XY returns the grid points for grid 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int GN, the grid size.
//
//    Output, double GX[GN], GY[GN], the grid coordinates.
//
{
  double x[33] = {
   0.05,  0.00,  0.00,  0.00,  0.10,
   0.10,  0.15,  0.20,  0.25,  0.30,
   0.35,  0.50,  0.50,  0.55,  0.60,
   0.60,  0.60,  0.65,  0.70,  0.70,
   0.70,  0.75,  0.75,  0.75,  0.80,
   0.80,  0.85,  0.90,  0.90,  0.95,
   1.00,  1.00,  1.00 };
  double y[33] = {
   0.45,  0.50,  1.00,  0.00,  0.15,
   0.75,  0.30,  0.10,  0.20,  0.35,
   0.85,  0.00,  1.00,  0.95,  0.25,
   0.65,  0.85,  0.70,  0.20,  0.65,
   0.90,  0.10,  0.35,  0.85,  0.40,
   0.65,  0.25,  0.35,  0.80,  0.90,
   0.00,  0.50,  1.00 };

  r8vec_copy ( gn, x, gx );
  r8vec_copy ( gn, y, gy );

  return;
}
//****************************************************************************80

int g03_size ( )

//****************************************************************************80
//
//  Purpose:
//
//    G03_SIZE returns the size for grid 3.
//
//  Discussion:
//
//    The "grid size" is simply the number of points in the grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GN, the grid size.
//
{
  int gn;

  gn = 25;

  return gn;
}
//****************************************************************************80

string g03_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    G03_TITLE returns the title for grid 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string GT, the grid title.
//
{
  string gt;

  gt = "Lawson's 25 node set";

  return gt;
}
//****************************************************************************80

void g03_xy ( int gn, double gx[], double gy[] )

//****************************************************************************80
//
//  Purpose:
//
//    G03_XY returns the grid points for grid 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int GN, the grid size.
//
//    Output, double GX[GN], GY[GN], the grid coordinates.
//
{
  double x[25] = {
   0.13750,   0.91250,   0.71250,   0.22500,  -0.05000,
   0.47500,   0.05000,   0.45000,   1.08750,   0.53750,
  -0.03750,   0.18750,   0.71250,   0.85000,   0.70000,
   0.27500,   0.45000,   0.81250,   0.45000,   1.00000,
   0.50000,   0.18750,   0.58750,   1.05000,   0.10000 };
  double y[25] = {
   0.97500,   0.98750,   0.76250,   0.83750,   0.41250,
   0.63750,  -0.05000,   1.03750,   0.55000,   0.80000,
   0.75000,   0.57500,   0.55000,   0.43750,   0.31250,
   0.42500,   0.28750,   0.18750,  -0.03750,   0.26250,
   0.46250,   0.26250,   0.12500,  -0.06125,   0.11250 };

  r8vec_copy ( gn, x, gx );
  r8vec_copy ( gn, y, gy );

  return;
}
//****************************************************************************80

int g04_size ( )

//****************************************************************************80
//
//  Purpose:
//
//    G04_SIZE returns the size for grid 4.
//
//  Discussion:
//
//    The "grid size" is simply the number of points in the grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GN, the grid size.
//
{
  int gn;

  gn = 100;

  return gn;
}
//****************************************************************************80

string g04_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    G04_TITLE returns the title for grid 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string GT, the grid title.
//
{
  string gt;

  gt = "Random 100 node set";

  return gt;
}
//****************************************************************************80

void g04_xy ( int gn, double gx[], double gy[] )

//****************************************************************************80
//
//  Purpose:
//
//    G04_XY returns the grid points for grid 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int GN, the grid size.
//
//    Output, double GX[GN], GY[GN], the grid coordinates.
//
{
  double x[100] = {
   0.0096326,  0.0216348,  0.0298360,  0.0417447,  0.0470462,
   0.0562965,  0.0646857,  0.0740377,  0.0873907,  0.0934832,
   0.1032216,  0.1110176,  0.1181193,  0.1251704,  0.1327330,
   0.1439536,  0.1564861,  0.1651043,  0.1786039,  0.1886405,
   0.2016706,  0.2099886,  0.2147003,  0.2204141,  0.2343715,
   0.2409660,  0.2527740,  0.2570839,  0.2733365,  0.2853833,
   0.2901755,  0.2964854,  0.3019725,  0.3125695,  0.3307163,
   0.3378504,  0.3439061,  0.3529922,  0.3635507,  0.3766172,
   0.3822429,  0.3869838,  0.3973137,  0.4170708,  0.4255588,
   0.4299218,  0.4372839,  0.4705033,  0.4736655,  0.4879299,
   0.4940260,  0.5055324,  0.5162593,  0.5219219,  0.5348529,
   0.5483213,  0.5569571,  0.5638611,  0.5784908,  0.5863950,
   0.5929148,  0.5987839,  0.6117561,  0.6252296,  0.6331381,
   0.6399048,  0.6488972,  0.6558537,  0.6677405,  0.6814074,
   0.6887812,  0.6940896,  0.7061687,  0.7160957,  0.7317445,
   0.7370798,  0.7462030,  0.7566957,  0.7699998,  0.7879347,
   0.7944014,  0.8164468,  0.8192794,  0.8368405,  0.8500993,
   0.8588255,  0.8646496,  0.8792329,  0.8837536,  0.8900077,
   0.8969894,  0.9044917,  0.9083947,  0.9203972,  0.9347906,
   0.9434519,  0.9490328,  0.9569571,  0.9772067,  0.9983493 };
  double y[100] = {
   0.3083158,  0.2450434,  0.8613847,  0.0977864,  0.3648355,
   0.7156339,  0.5311312,  0.9755672,  0.1781117,  0.5452797,
   0.1603881,  0.7837139,  0.9982015,  0.6910589,  0.1049580,
   0.8184662,  0.7086405,  0.4456593,  0.1178342,  0.3189021,
   0.9668446,  0.7571834,  0.2016598,  0.3232444,  0.4368583,
   0.8907869,  0.0647260,  0.5692618,  0.2947027,  0.4332426,
   0.3347464,  0.7436284,  0.1066265,  0.8845357,  0.5158730,
   0.9425637,  0.4799701,  0.1783069,  0.1146760,  0.8225797,
   0.2270688,  0.4073598,  0.8875080,  0.7631616,  0.9972804,
   0.4959884,  0.3410421,  0.2498120,  0.6409007,  0.1058690,
   0.5411969,  0.0089792,  0.8784268,  0.5515874,  0.4038952,
   0.1654023,  0.2965158,  0.3660356,  0.0366554,  0.9502420,
   0.2638101,  0.9277386,  0.5377694,  0.7374676,  0.4674627,
   0.9186109,  0.0416884,  0.1291029,  0.6763676,  0.8444238,
   0.3273328,  0.1893879,  0.0645923,  0.0180147,  0.8904992,
   0.4160648,  0.4688995,  0.2174508,  0.5734231,  0.8853319,
   0.8018436,  0.6388941,  0.8931002,  0.1000558,  0.2789506,
   0.9082948,  0.3259159,  0.8318747,  0.0508513,  0.9708450,
   0.5120548,  0.2859716,  0.9581641,  0.6183429,  0.3779934,
   0.4010423,  0.9478657,  0.7425486,  0.8883287,  0.5496750 };

  r8vec_copy ( gn, x, gx );
  r8vec_copy ( gn, y, gy );

  return;
}
//****************************************************************************80

int g05_size ( )

//****************************************************************************80
//
//  Purpose:
//
//    G05_SIZE returns the size for grid 5.
//
//  Discussion:
//
//    The "grid size" is simply the number of points in the grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GN, the grid size.
//
{
  int gn;

  gn = 81;

  return gn;
}
//****************************************************************************80

string g05_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    G05_TITLE returns the title for grid 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string GT, the grid title.
//
{
  string gt;

  gt = "Gridded 81 node set";

  return gt;
}
//****************************************************************************80

void g05_xy ( int gn, double gx[], double gy[] )

//****************************************************************************80
//
//  Purpose:
//
//    G05_XY returns the grid points for grid 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int GN, the grid size.
//
//    Output, double GX[GN], GY[GN], the grid coordinates.
//
{
  double x[81] = {
   0.000,  0.000,  0.000,  0.000,  0.000,
   0.000,  0.000,  0.000,  0.000,  0.125,
   0.125,  0.125,  0.125,  0.125,  0.125,
   0.125,  0.125,  0.125,  0.250,  0.250,
   0.250,  0.250,  0.250,  0.250,  0.250,
   0.250,  0.250,  0.375,  0.375,  0.375,
   0.375,  0.375,  0.375,  0.375,  0.375,
   0.375,  0.500,  0.500,  0.500,  0.500,
   0.500,  0.500,  0.500,  0.500,  0.500,
   0.625,  0.625,  0.625,  0.625,  0.625,
   0.625,  0.625,  0.625,  0.625,  0.750,
   0.750,  0.750,  0.750,  0.750,  0.750,
   0.750,  0.750,  0.750,  0.875,  0.875,
   0.875,  0.875,  0.875,  0.875,  0.875,
   0.875,  0.875,  1.000,  1.000,  1.000,
   1.000,  1.000,  1.000,  1.000,  1.000,
   1.000 };
  double y[81] = {
   0.000,  0.125,  0.250,  0.375,  0.500,
   0.625,  0.750,  0.875,  1.000,  0.000,
   0.125,  0.250,  0.375,  0.500,  0.625,
   0.750,  0.875,  1.000,  0.000,  0.125,
   0.250,  0.375,  0.500,  0.625,  0.750,
   0.875,  1.000,  0.000,  0.125,  0.250,
   0.375,  0.500,  0.625,  0.750,  0.875,
   1.000,  0.000,  0.125,  0.250,  0.375,
   0.500,  0.625,  0.750,  0.875,  1.000,
   0.000,  0.125,  0.250,  0.375,  0.500,
   0.625,  0.750,  0.875,  1.000,  0.000,
   0.125,  0.250,  0.375,  0.500,  0.625,
   0.750,  0.875,  1.000,  0.000,  0.125,
   0.250,  0.375,  0.500,  0.625,  0.750,
   0.875,  1.000,  0.000,  0.125,  0.250,
   0.375,  0.500,  0.625,  0.750,  0.875,
   1.000 };

  r8vec_copy ( gn, x, gx );
  r8vec_copy ( gn, y, gy );

  return;
}
