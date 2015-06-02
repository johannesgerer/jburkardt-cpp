# include <cstdlib>
# include <iostream>
# include <cmath>

using namespace std;

int main ( );
void sor_test01 ( double w );

# include "gnuplot_i.hpp"
# include "sor.hpp"

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SOR_PRB.
//
//  Discussion:
//
//    SOR_PRB tests the SOR library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  double w;

  timestamp ( );
  cout << "\n";
  cout << "SOR_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the SOR library.\n";

  w = 0.5;
  sor_test01 ( w );
  sleep ( 2 );

  w = 1.0;
  sor_test01 ( w );
  sleep ( 2 );

  w = 1.5;
  sor_test01 ( w );
  sleep ( 2 );
//
//  Terminate.
//
  cout << "\n";
  cout << "SOR_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void sor_test01 ( double w )

//****************************************************************************80
//
//  Purpose:
//
//    SOR_TEST01 tests SOR1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double W, the relaxation factor.
//    0 < W < 2 is required.
//
{
  double *a;
  double *b;
  int i;
  int it;
  int it_num;
  int j;
  double *m_plot;
  double *ml_plot;
  int n;
  double *r_plot;
  double *rl_plot;
  double *s_plot;
  double t;
  gnuplot_ctrl *window1;
  gnuplot_ctrl *window2;
  double *x;
  double *x_exact;
  double *x_new;
  double *x_plot;

  cout << "\n";
  cout << "SOR1_TEST01:\n";
  cout << "  Relaxation parameter W = " << w << "\n";

  it_num = 2000;
  n = 33;
//
//  Set the matrix A.
//
  a = dif2 ( n, n );
//
//  Determine the right hand side vector B.
//
  x_exact = new double[n];
  for ( i = 0; i < n; i++ )
  {
    t = ( double ) i / ( double ) ( n - 1 );
    x_exact[i] = exp ( t ) * ( t - 1 ) * t;
//   x_exact[i] = ( double ) ( i + 1 );
  }
  b = r8mat_mv ( n, n, a, x_exact );
//
//  Set the initial estimate for the solution.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }
//
//  Allocate plot arrays.
//
  m_plot = new double[it_num+1];
  r_plot = new double[it_num+1];
  s_plot = new double[it_num+1];
  x_plot = new double[n*(it_num+1)];
//
//  Initialize plot arrays.
//
  r_plot[0] = r8mat_residual_norm ( n, n, a, x, b );
  m_plot[0] = 1.0;
  for ( i = 0; i < n; i++ )
  {
    x_plot[i+0*n] = x[i];
  }
  for ( j = 0; j <= it_num; j++ )
  {
    s_plot[j] = ( double ) j;
  }
//
//  Carry out the iteration.
//
  for ( it = 1; it <= it_num; it++ )
  {
    x_new = sor1 ( n, a, b, x, w );

    r_plot[it] = r8mat_residual_norm ( n, n, a, x_new, b );
//
//  Compute the average point motion.
//
    m_plot[it] = r8vec_diff_norm_squared ( n, x, x_new ) / ( double ) n;
//
//  Update the solution
//
    r8vec_copy ( n, x_new, x );

    for ( i = 0; i < n; i++ )
    {
      x_plot[i+0*n] = x[i];
    }

    delete [] x_new;
  }
    r8vec_print ( n, x, "Solution" );
//
//  Plot the residual.
//
  rl_plot = new double[it_num+1];
  for ( j = 0; j <= it_num; j++ )
  {
    rl_plot[j] = log ( r_plot[j] );
  }
  window1 = gnuplot_init ( );
  gnuplot_setstyle ( window1, "lines" );
  gnuplot_cmd ( window1, "set grid" );
  gnuplot_set_xlabel ( window1, "Step" );
  gnuplot_set_ylabel ( window1, "Log ( Residual )" );
  gnuplot_setstyle ( window1, "lines" );
  gnuplot_plot1d_var2v ( window1, s_plot, rl_plot, it_num + 1, "Log(Residual)" );
  sleep ( 1 );
//
//  Plot the average point motion.
//
  ml_plot = new double[it_num+1];
  for ( j = 0; j <= it_num; j++ )
  {
    ml_plot[j] = log ( m_plot[j] );
//  cout << j << "  " << ml_plot[j] << "\n";
  }
  window2 = gnuplot_init ( );
  gnuplot_setstyle ( window2, "lines" );
  gnuplot_cmd ( window2, "set grid" );
  gnuplot_set_xlabel ( window2, "Step" );
  gnuplot_set_ylabel ( window2, "Log ( Motion )" );
  gnuplot_setstyle ( window2, "lines" );
  gnuplot_plot1d_var2v ( window2, s_plot, ml_plot, it_num + 1, "Log(Motion)" );
  sleep ( 1 );
//
//  Plot the evolution of the locations of the generators.
//
//figure ( 3 )

//y = ( 0 : it_num );
//for k = 1 : n
//  plot ( x_plot(k,1:it_num+1), y )
//  hold on;
//end
//grid on
//hold off;

//title ( "Generator evolution." );
//xlabel ( "Generator positions" );
//ylabel ( "Iterations" ); 

  gnuplot_close ( window1 );
  gnuplot_close ( window2 );

  delete [] a;
  delete [] b;
  delete [] m_plot;
  delete [] ml_plot;
  delete [] r_plot;
  delete [] rl_plot;
  delete [] s_plot;
  delete [] x;
  delete [] x_exact;
  delete [] x_plot;

  return;
}
