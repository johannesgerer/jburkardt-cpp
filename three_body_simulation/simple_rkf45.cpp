# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstring>
# include <cmath>
# include <ctime>

using namespace std;

# include "rkf45.hpp"

int main ( );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void simple_rkf45_run ( );
void simple_f ( double t, double y[], double yp[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SIMPLE_RKF45.
//
//  Discussion:
//
//    SIMPLE_RKF45 uses RKF45 as an integrator for the simple version
//    of the three-body problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "SIMPLE_RKF45\n";
  cout << "  C++ version\n";
  cout << "  Simulate the behavior of three bodies which are\n";
  cout << "  constrained to lie in a plane, moving under the\n";
  cout << "  influence of gravity.\n";
  cout << " \n";
  cout << "  Use RKF45 for the ODE integrator.\n";

  simple_rkf45_run ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SIMPLE_RKF45\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
//****************************************************************************80

void simple_rkf45_run ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLE_RKF45_RUN runs the simple three body ODE system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  double abserr;
  int flag;
  int i;
  double m0;
  double m1;
  double m2;
  int neqn = 12;
  double relerr;
  int step;
  int step_num = 630;
  double t;
  string t_filename = "simple_rkf45_t.txt";
  double t_out;
  double t_start;
  double t_stop;
  double *ts;
  double *y;
  string y_filename = "simple_rkf45_y.txt";
  double *yp;
  double *ys;

  ts = new double[step_num+1];
  y = new double[neqn];
  yp = new double[neqn];
  ys = new double[neqn*(step_num+1)];

  cout << "\n";
  cout << "SIMPLE_RKF45_RUN\n";
  cout << "  Simulate the planar three-body problem as an ODE system\n";
  cout << "  using RKF45 for the ODE integration.\n";

  m0 = 5.0;
  m1 = 3.0;
  m2 = 4.0;

  abserr = 1.0E-10;
  relerr = 1.0E-10;

  flag = 1;

  t_start = 0.0;
  t_stop = 63.0;

  t = 0.0;
  t_out = 0.0;

  y[0] =  1.0;
  y[1] = -1.0;
  y[2] =  0.0;
  y[3] =  0.0;
  y[4] =  1.0;
  y[5] =  3.0;
  y[6] =  0.0;
  y[7] =  0.0;
  y[8] = -2.0;
  y[9] = -1.0;
  y[10] = 0.0;
  y[11] = 0.0;

  simple_f ( t, y, yp );

  for ( i = 0; i < neqn; i++ )
  {
    ys[i+0*neqn] = y[i];
  }
  ts[0] = t;

  for ( step = 1; step <= step_num; step++ )
  {
    t = ( ( double ) ( step_num - step + 1 ) * t_start 
        + ( double ) (            step - 1 ) * t_stop ) 
        / ( double ) ( step_num            );

    t_out = ( ( double ) ( step_num - step ) * t_start 
            + ( double ) (            step ) * t_stop ) 
            / ( double ) ( step_num        );

    flag = r8_rkf45 ( simple_f, neqn, y, yp, &t, t_out, &relerr, abserr, flag );

    if ( abs ( flag ) != 2 )
    {
      cout << "\n";
      cout << "SIMPLE_RKF45_RUN - Warning\n";
      cout << "  Output value of FLAG = " << flag 
           << " at output time T_OUT = " << t_out << "\n";
    }

    if ( flag == 7 )
    {
      flag = 2;
    }

    for ( i = 0; i < neqn; i++ )
    {
      ys[i+step*neqn] = y[i];
    }
    ts[step] = t_out;
  }

  r8mat_write ( t_filename, 1, step_num + 1, ts );
  r8mat_write ( y_filename, neqn, step_num + 1, ys );

  cout << "\n";
  cout << "SIMPLE_RKF45_RUN:\n";
  cout << "  Time data written to \"" << t_filename << "\".\n";
  cout << "  Solution data written to \"" << y_filename << "\".\n";

  delete [] ts;
  delete [] y;
  delete [] yp;
  delete [] ys;

  return;
}
//****************************************************************************80

void simple_f ( double t, double y[], double yp[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLE_F returns the right hand side of the three body ODE system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2012
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
//    Output, double YP[NEQN], the value of the derivatives.
//
{
  double m0;
  double m1;
  double m2;
  double n0;
  double n1;
  double n2;
  int neqn = 12;
  double x0;
  double x1;
  double x2;
  double y0;
  double y1;
  double y2;

  m0 = 5.0;
  m1 = 3.0;
  m2 = 4.0;

  x0 = y[0];
  y0 = y[1];

  x1 = y[4];
  y1 = y[5];

  x2 = y[8];
  y2 = y[9];

  n0 = sqrt ( pow ( pow ( x2 - x1, 2 ) + pow ( y2 - y1, 2 ), 3 ) );
  n1 = sqrt ( pow ( pow ( x0 - x2, 2 ) + pow ( y0 - y2, 2 ), 3 ) );
  n2 = sqrt ( pow ( pow ( x1 - x0, 2 ) + pow ( y1 - y0, 2 ), 3 ) ); 

  yp[0]  =  y[2];
  yp[1]  =  y[3];
  yp[2]  = - m1 * ( x0 - x1 ) / n2 - m2 * ( x0 - x2 ) / n1;
  yp[3]  = - m1 * ( y0 - y1 ) / n2 - m2 * ( y0 - y2 ) / n1;
  yp[4]  =  y[6];
  yp[5]  =  y[7];
  yp[6]  = - m2 * ( x1 - x0 ) / n0 - m0 * ( x1 - x2 ) / n2;
  yp[7]  = - m2 * ( y1 - y0 ) / n0 - m0 * ( y1 - y2 ) / n2;
  yp[8]  = y[10];
  yp[9]  = y[11];
  yp[10] = - m0 * ( x2 - x0 ) / n1 - m1 * ( x2 - x1 ) / n0;
  yp[11] = - m0 * ( y2 - y0 ) / n1 - m1 * ( y2 - y1 ) / n0;

  return;
}
