# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "interp.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( int data_num );
void test04 ( int data_num );
double *f_runge ( int m, int n, double x[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for INTERP_PRB.
//
//  Discussion:
//
//    INTERP_PRB tests the INTERP library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
{
  int data_num;

  timestamp ( );
  cout << "\n";
  cout << "INTERP_PRB\n";
  cout << "  C++ version:\n";
  cout << "  Test the INTERP library.\n";

  test01 ( );

  test02 ( );

  data_num = 6;
  test03 ( data_num );

  data_num = 11;
  test03 ( data_num );

  data_num = 6;
  test04 ( data_num );

  data_num = 11;
  test04 ( data_num );
//
//  Terminate.
//
  cout << "\n";
  cout << "INTERP_PRB\n";
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
//    TEST01 tests INTERP_NEAREST on 1-dimensional data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
{
  int after;
  int before;
  int data_num = 11;
  int fat;
  int i;
  int interp;
  int interp_num;
  int j;
  int m = 1;
  double p;
  double *p_data;
  double *p_interp;
  double *p_value;
  double t;
  double *t_data;
  double *t_interp;
  double t_max;
  double t_min;

  cout << " \n";
  cout << "TEST01\n";
  cout << "  INTERP_NEAREST evaluates a nearest-neighbor interpolant.\n";
  cout << " \n";
  cout << "  In this example, the function we are interpolating is\n";
  cout << "  Runge's function, with Chebyshev knots.\n";

  t_min = -1.0;
  t_max = +1.0;

  t_data = cc_abscissas_ab ( t_min, t_max, data_num );

  p_data = f_runge ( m, data_num, t_data );

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Spatial dimension =     " << m << "\n";
  cout << "  Number of data values = " << data_num << "\n";
  cout << "\n";
  cout << "       T_data        P_data\n";
  cout << "\n";
  for ( j = 0; j < data_num; j++ )
  {
    cout << "  " << setw(14) << t_data[j]
         << "  " << setw(14) << p_data[0+j*m] << "\n";
  }
//
//  Our interpolation values will include the original T values, plus
//  3 new values in between each pair of original values.
//
  before = 4;
  fat = 3;
  after = 2;

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after;

  t_interp = r8vec_expand_linear2 ( data_num, t_data, before, fat, after );

  p_interp = interp_nearest ( m, data_num, t_data, p_data, interp_num, t_interp );

  p_value = f_runge ( m, interp_num, t_interp );

  cout << "\n";
  cout << "  Interpolation:\n";
  cout << "\n";
  cout << "    T_interp      P_interp        P_exact        Error\n";
  cout << "\n";

  for ( j = 0; j < interp_num; j++ )
  {
    cout << "  " << setw(10) << t_interp[j]
         << "  " << setw(14) << p_interp[0+j*m]
         << "  " << setw(14) << p_value[0+j*m]
         << "  " << setw(10) << p_interp[0+j*m] - p_value[0+j*m] << "\n";
  }

  delete [] p_data;
  delete [] p_interp;
  delete [] p_value;
  delete [] t_data;
  delete [] t_interp;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests INTERP_LINEAR on 1-dimensional data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
{
  int after;
  int before;
  int data_num = 11;
  int fat;
  int i;
  int interp;
  int interp_num;
  int j;
  int m = 1;
  double p;
  double *p_data;
  double *p_interp;
  double *p_value;
  double t;
  double *t_data;
  double *t_interp;
  double t_max;
  double t_min;

  cout << " \n";
  cout << "TEST02\n";
  cout << "  INTERP_LINEAR evaluates a piecewise linear spline.\n";
  cout << " \n";
  cout << "  In this example, the function we are interpolating is\n";
  cout << "  Runge's function, with evenly spaced knots.\n";

  t_min = -1.0;
  t_max = +1.0;

  t_data = ncc_abscissas_ab ( t_min, t_max, data_num );

  p_data = f_runge ( m, data_num, t_data );

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Spatial dimension =     " << m << "\n";
  cout << "  Number of data values = " << data_num << "\n";
  cout << "\n";
  cout << "       T_data        P_data\n";
  cout << "\n";
  for ( j = 0; j < data_num; j++ )
  {
    cout << "  " << setw(14) << t_data[j]
         << "  " << setw(14) << p_data[0+j*m] << "\n";
  }
//
//  Our interpolation values will include the original T values, plus
//  3 new values in between each pair of original values.
//
  before = 4;
  fat = 3;
  after = 2;

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after;

  t_interp = r8vec_expand_linear2 ( data_num, t_data, before, fat, after );

  p_interp = interp_linear ( m, data_num, t_data, p_data, interp_num, 
    t_interp );

  p_value = f_runge ( m, interp_num, t_interp );

  cout << "\n";
  cout << "  Interpolation:\n";
  cout << "\n";
  cout << "    T_interp      P_interp        P_exact        Error\n";
  cout << "\n";

  for ( j = 0; j < interp_num; j++ )
  {
    cout << "  " << setw(10) << t_interp[j]
         << "  " << setw(14) << p_interp[0+j*m]
         << "  " << setw(14) << p_value[0+j*m]
         << "  " << setw(10) << p_interp[0+j*m] - p_value[0+j*m] << "\n";
  }

  delete [] p_data;
  delete [] p_interp;
  delete [] p_value;
  delete [] t_data;
  delete [] t_interp;

  return;
}
//****************************************************************************80

void test03 ( int data_num )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests INTERP_LAGRANGE on 1-dimensional data, equally spaced data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the number of data values.
//
{
  int after;
  int before;
  int fat;
  int i;
  int interp;
  int interp_num;
  int j;
  int m = 1;
  double p;
  double *p_data;
  double *p_interp;
  double *p_value;
  double t;
  double *t_data;
  double *t_interp;
  double t_max;
  double t_min;

  cout << " \n";
  cout << "TEST03\n";
  cout << "  INTERP_LAGRANGE evaluates a polynomial interpolant.\n";
  cout << "  In this example, the function we are interpolating is\n";
  cout << "  Runge's function, with evenly spaced knots.\n";

  t_min = -1.0;
  t_max = +1.0;

  t_data = ncc_abscissas_ab ( t_min, t_max, data_num );

  p_data = f_runge ( m, data_num, t_data );

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Spatial dimension =     " << m << "\n";
  cout << "  Number of data values = " << data_num << "\n";
  cout << "\n";
  cout << "       T_data        P_data\n";
  cout << "\n";
  for ( j = 0; j < data_num; j++ )
  {
    cout << "  " << setw(14) << t_data[j]
         << "  " << setw(14) << p_data[0+j*m] << "\n";
  }
//
//  Our interpolation values will include the original T values, plus
//  3 new values in between each pair of original values.
//
  before = 4;
  fat = 3;
  after = 2;

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after;

  t_interp = r8vec_expand_linear2 ( data_num, t_data, before, fat, after );

  p_interp = interp_lagrange ( m, data_num, t_data, p_data, interp_num, 
    t_interp );

  p_value = f_runge ( m, interp_num, t_interp );

  cout << "\n";
  cout << "  Interpolation:\n";
  cout << "\n";
  cout << "    T_interp      P_interp        P_exact        Error\n";
  cout << "\n";

  for ( j = 0; j < interp_num; j++ )
  {
    cout << "  " << setw(10) << t_interp[j]
         << "  " << setw(14) << p_interp[0+j*m]
         << "  " << setw(14) << p_value[0+j*m]
         << "  " << setw(10) << p_interp[0+j*m] - p_value[0+j*m] << "\n";
  }

  delete [] p_data;
  delete [] p_interp;
  delete [] p_value;
  delete [] t_data;
  delete [] t_interp;

  return;
}
//****************************************************************************80

void test04 ( int data_num )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests INTERP_LAGRANGE on 1-dimensional data, Clenshaw-Curtis data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the number of data values.
//
{
  int after;
  int before;
  int fat;
  int i;
  int interp;
  int interp_num;
  int j;
  int m = 1;
  double p;
  double *p_data;
  double *p_interp;
  double *p_value;
  double t;
  double *t_data;
  double *t_interp;
  double t_max;
  double t_min;

  cout << " \n";
  cout << "TEST04\n";
  cout << "  INTERP_LAGRANGE evaluates a polynomial interpolant.\n";
  cout << " \n";
  cout << "  In this example, the function we are interpolating is\n";
  cout << "  Runge's function, with Clenshaw Curtis knots.\n";

  t_min = -1.0;
  t_max = +1.0;

  t_data = cc_abscissas_ab ( t_min, t_max, data_num );

  p_data = f_runge ( m, data_num, t_data );

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Spatial dimension =     " << m << "\n";
  cout << "  Number of data values = " << data_num << "\n";
  cout << "\n";
  cout << "       T_data        P_data\n";
  cout << "\n";
  for ( j = 0; j < data_num; j++ )
  {
    cout << "  " << setw(14) << t_data[j]
         << "  " << setw(14) << p_data[0+j*m] << "\n";
  }
//
//  Our interpolation values will include the original T values, plus
//  3 new values in between each pair of original values.
//
  before = 4;
  fat = 3;
  after = 2;

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after;

  t_interp = r8vec_expand_linear2 ( data_num, t_data, before, fat, after );

  p_interp = interp_lagrange ( m, data_num, t_data, p_data, interp_num, 
    t_interp );

  p_value = f_runge ( m, interp_num, t_interp );

  cout << "\n";
  cout << "  Interpolation:\n";
  cout << "\n";
  cout << "    T_interp      P_interp        P_exact        Error\n";
  cout << "\n";

  for ( j = 0; j < interp_num; j++ )
  {
    cout << "  " << setw(10) << t_interp[j]
         << "  " << setw(14) << p_interp[0+j*m]
         << "  " << setw(14) << p_value[0+j*m]
         << "  " << setw(10) << p_interp[0+j*m] - p_value[0+j*m] << "\n";
  }

  delete [] p_data;
  delete [] p_interp;
  delete [] p_value;
  delete [] t_data;
  delete [] t_interp;

  return;
}
//****************************************************************************80

double *f_runge ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    F_RUNGE evaluates the Runge function.
//
//  Discussion:
//
//    Interpolation of the Runge function at evenly spaced points in [-1,1]
//    is a common test.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double F_RUNGE[N], the function values.
//
{
  double *f;
  int i;
  int j;
  double t;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    t = 0.0;
    for ( i = 0; i < m; i++ )
    {
      t = t + pow ( x[i+j*m], 2 );
    }
    f[j] = 1.0 / ( 1.0 + 25.0 * t );
  }

  return f;
}
