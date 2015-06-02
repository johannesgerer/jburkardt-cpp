# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "svd_snowfall.hpp"

int main ( );
void svd_snowfall_test02 ( int m, int n, double x[] );
void svd_snowfall_test03 ( int m, int n, double x[] );
void svd_snowfall_test04 ( int m, int n, double x[] );
void svd_snowfall_test05 ( int m, int n, double x[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SVD_SNOWFALL_PRB.
//
//  Discussion:
//
//    SVD_SNOWFALL_PRB tests the SVD_SNOWFALL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 May 2013
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "snowfall.txt";
  int m;
  int n;
  double *x;

  timestamp ( );
  cout << "\n";
  cout << "SVD_SNOWFALL_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SVD_SNOWFALL library.\n";
//
//  Retrieve the data.
//  It's really easier to do this in the main program.
//
  cout << "\n";
  cout << "SVD_SNOWFALL_TEST01\n";
  cout << "  Read, process, and return snowfall data in \"" << filename << "\".\n";
//
//  Determine the size of the data.
//
  r8mat_header_read ( filename, &m, &n );

  cout << "\n";
  cout << "  Number of data rows    M = " << m << "\n";
  cout << "  Number of data columns N = " << n << "\n";

  x = r8mat_data_read ( filename, m, n );

  cout << "\n";
  cout << "  Data has been read from the file.\n";

  svd_snowfall_test02 ( m, n, x );
  svd_snowfall_test03 ( m, n, x );
  svd_snowfall_test04 ( m, n, x );
  svd_snowfall_test05 ( m, n, x );
//
//  Free memory.
//
  free ( x );
//
//  Terminate.
//
  cout << "\n";
  cout << "SVD_SNOWFALL_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void svd_snowfall_test02 ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SVD_SNOWFALL_TEST02 looks at the singular values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double X[M*N], the snowfall data.
//
{
  string command_filename = "singular_values_commands.txt";
  ofstream command;
  string data_filename = "singular_values_data.txt";
  ofstream data;
  double *e;
  double *e_cum;
  double e_sum;
  int i;
  int mn;
  double *s;
  double *s_diag;
  double *u;
  double *v;

  cout << "\n";
  cout << "SVD_SNOWFALL_TEST02\n";
  cout << "  Look at the singular values.\n";
  cout << "  If the singular values are close, then the data is\n";
  cout << "  well spread out.  If the singular values decay rapidly,\n";
  cout << "  then the data exhibits patterns, or is constrained to\n";
  cout << "  a lower-dimensional subspace.\n";
//
//  Compute the SVD.
//
  u = new double[m*m];
  s = new double[m*n];
  v = new double[n*n];

  r8mat_svd_linpack ( m, n, x, u, s, v );
//
//  Extract the diagonal of S.
//
  mn = i4_min ( m, n );
  s_diag = ( double * ) malloc ( mn * sizeof ( double ) );

  for ( i = 0; i < mn; i++ )
  {
    s_diag[i] = s[i+i*m];
  }
//
//  Print the singular values.
//
  r8vec_print ( mn, s_diag, "  The singular values:" );
//
//  Plot the singular values.
//
  data.open ( data_filename.c_str ( ) );
  for ( i = 0; i < mn; i++ )
  {
    data << "  " << setw(4) << i
         << "  " << setw(14) << s_diag[i] << "\n";
  }
  data.close ( );
  cout << "\n";
  cout << "  Created data file \"" << data_filename << "\".\n";

  command.open ( command_filename.c_str ( ) );
  command << "# " << command_filename << "\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < " << command_filename << "\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'singular_values.png'\n";
  command << "set xlabel 'Index I'\n";
  command << "set ylabel 'S(I)'\n";
  command << "set title 'Snowfall Singular Values'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'blue'\n"; 
  command << "quit\n";
  command.close ( );
  cout << "  Created command file \"" << command_filename << "\".\n";
//
//  Print the cumulative "energy" of the singular values.
//
  e = new double[mn];

  for ( i = 0; i < mn; i++ )
  {
    e[i] = pow ( s_diag[i], 2 );
  }
  e_sum = r8vec_sum ( mn, e );
  for ( i = 0; i < mn; i++ )
  {
    e[i] = e[i] / e_sum;
  }
  e_cum = r8vec_cum0_new ( mn, e );

  r8vec_print ( mn + 1, e_cum, "  The cumulative energy:" );

  delete [] e;
  delete [] e_cum;
  delete [] s;
  delete [] s_diag;
  delete [] u;
  delete [] v;

  return;
}
//****************************************************************************80

void svd_snowfall_test03 ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SVD_SNOWFALL_TEST03 computes low rank approximations to the matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double X[M*N], the snowfall data.
//
{
  double *a1;
  double *a2;
  double *a3;
  double *a4;
  double *a5;
  string command_filename = "approx_commands.txt";
  ofstream command;
  string data_filename = "approx_data.txt";
  ofstream data;
  int i;
  double *s;
  double *u;
  double *v;

  cout << "\n";
  cout << "SVD_SNOWFALL_TEST03\n";
  cout << "  Compute the rank 1 through rank 5 approximations to the data.\n";
  cout << "  Compare each of these to the 2012 snowfall data.\n";
//
//  Compute the SVD.
//
  u = new double[m*m];
  s = new double[m*n];
  v = new double[n*n];

  r8mat_svd_linpack ( m, n, x, u, s, v );
//
//  Form the rank 1, 2, 3, 4, 5 approximants to A.
//
  a1 = r8mat_svd_low_rank ( m, n, 1, u, s, v );
  a2 = r8mat_svd_low_rank ( m, n, 2, u, s, v );
  a3 = r8mat_svd_low_rank ( m, n, 3, u, s, v );
  a4 = r8mat_svd_low_rank ( m, n, 4, u, s, v );
  a5 = r8mat_svd_low_rank ( m, n, 5, u, s, v );
//
//  Column 1 of X is the 2012 snowfall.
//  Column 1 of A1 is the rank 1 approximant to 2012 snowfall.
//
  data.open ( data_filename.c_str ( ) );
  for ( i = 0; i < m; i++ )
  {
    data << "  " << setw(4) << i
         << "  " << setw(14) << x[i+0*m]
         << "  " << setw(14) << a1[i+0*m]
         << "  " << setw(14) << a2[i+0*m]
         << "  " << setw(14) << a3[i+0*m]
         << "  " << setw(14) << a4[i+0*m]
         << "  " << setw(14) << a5[i+0*m] << "\n";
  }
  data.close ( );
  cout << "  Created data file \"" << data_filename << "\".\n";

  command.open ( command_filename.c_str ( ) );
  command << "# " << command_filename << "\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < " << command_filename << "\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'approx0.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title '2012 Snowfall'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'blue'\n";

  command << "set output 'approx1.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Rank 1 Approx to 2012 Snowfall'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'blue',\\\n";
  command << "     '" << data_filename << "' using 1:3 lw 3 linecolor rgb 'red'\n";

  command << "set output 'approx2.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Rank 2 Approx to 2012 Snowfall'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'blue',\\\n";
  command << "     '" << data_filename << "' using 1:3 lw 3 linecolor rgb 'gray',\\\n";
  command << "     '" << data_filename << "' using 1:4 lw 3 linecolor rgb 'red'\n";

  command << "set output 'approx3.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Rank 3 Approx to 2012 Snowfall'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'blue',\\\n";
  command << "     '" << data_filename << "' using 1:3 lw 3 linecolor rgb 'gray',\\\n";
  command << "     '" << data_filename << "' using 1:4 lw 3 linecolor rgb 'gray',\\\n";
  command << "     '" << data_filename << "' using 1:5 lw 3 linecolor rgb 'red'\n";

  command << "set output 'approx4.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Rank 4 Approx to 2012 Snowfall'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'blue',\\\n";
  command << "     '" << data_filename << "' using 1:3 lw 3 linecolor rgb 'gray',\\\n";
  command << "     '" << data_filename << "' using 1:4 lw 3 linecolor rgb 'gray',\\\n";
  command << "     '" << data_filename << "' using 1:5 lw 3 linecolor rgb 'gray',\\\n";
  command << "     '" << data_filename << "' using 1:6 lw 3 linecolor rgb 'red'\n";

  command << "set output 'approx5.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Rank 5 Approx to 2012 Snowfall'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'blue',\\\n";
  command << "     '" << data_filename << "' using 1:3 lw 3 linecolor rgb 'gray',\\\n";
  command << "     '" << data_filename << "' using 1:4 lw 3 linecolor rgb 'gray',\\\n";
  command << "     '" << data_filename << "' using 1:5 lw 3 linecolor rgb 'gray',\\\n";
  command << "     '" << data_filename << "' using 1:6 lw 3 linecolor rgb 'gray',\\\n";
  command << "     '" << data_filename << "' using 1:7 lw 3 linecolor rgb 'red'\n";

  command << "quit\n";
  command.close ( );
  cout << "  Created command file '" << command_filename << "'.\n";

  delete [] a1;
  delete [] a2;
  delete [] a3;
  delete [] a4;
  delete [] a5;
  delete [] s;
  delete [] u;
  delete [] v;

  return;
}
//****************************************************************************80

void svd_snowfall_test04 ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SVD_SNOWFALL_TEST04 looks at the first 6 modes in the U matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double X[M*N], the snowfall data.
//
{
  string command_filename = "umode_commands.txt";
  ofstream command;
  string data_filename = "umode_data.txt";
  ofstream data;
  int i;
  double *s;
  double *u;
  double *v;

  cout << "\n";
  cout << "SVD_SNOWFALL_TEST04\n";
  cout << "  Look at the first 6 modes in the U matrix.\n";
  cout << "  Each of these represents a pattern for snowfall over a year.\n";
  cout << "  The first mode is the pattern that is strongest in the data.\n";
//
//  Compute the SVD.
//
  u = new double[m*m];
  s = new double[m*n];
  v = new double[n*n];

  r8mat_svd_linpack ( m, n, x, u, s, v );
//
//  Normalize the patterns so that each column has maximum entry 1.
//
  r8col_normalize_li ( m, m, u );
//
//  Plot the U modes.
//
  data.open ( data_filename.c_str ( ) );
  for ( i = 0; i < m; i++ )
  {
    data << "  " << setw(4) << i
         << "  " << setw(14) << u[i+0*m]
         << "  " << setw(14) << u[i+1*m]
         << "  " << setw(14) << u[i+2*m]
         << "  " << setw(14) << u[i+3*m]
         << "  " << setw(14) << u[i+4*m]
         << "  " << setw(14) << u[i+5*m] << "\n";
  }
  data.close ( );
  cout << "  Created data file \"" << data_filename << "\".\n";

  command.open ( command_filename.c_str ( ) );
  command << "# " << command_filename << "\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < " << command_filename << "\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'umode1.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Monthly Snowfall Mode 1'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'blue'\n";

  command << "set output 'umode2.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Monthly Snowfall Mode 2'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:3 lw 3 linecolor rgb 'blue'\n";

  command << "set output 'umode3.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Monthly Snowfall Mode 3'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:4 lw 3 linecolor rgb 'blue'\n";
  command << "set output 'umode4.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Monthly Snowfall Mode 4'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:5 lw 3 linecolor rgb 'blue'\n";

  command << "set output 'umode5.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Monthly Snowfall Mode 5'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:6 lw 3 linecolor rgb 'blue'\n";

  command << "set output 'umode6.png'\n";
  command << "set xlabel 'Month'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Monthly Snowfall Mode 6'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:7 lw 3 linecolor rgb 'blue'\n";

  command << "quit\n";
  command.close ( );
  cout << "  Created command file '" << command_filename << "'.\n";

  delete [] s;
  delete [] u;
  delete [] v;

  return;
}
//****************************************************************************80

void svd_snowfall_test05 ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SVD_SNOWFALL_TEST05 looks at the first 6 modes in the V matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double X[M*N], the snowfall data.
//
{
  string command_filename = "vmode_commands.txt";
  ofstream command;
  string data_filename = "vmode_data.txt";
  ofstream data;
  int i;
  double *s;
  double *u;
  double *v;

  cout << "\n";
  cout << "SVD_SNOWFALL_TEST05\n";
  cout << "  Look at the first 6 modes in the V matrix.\n";
  cout << "  Each of these represents a pattern shared by all the months,\n";
  cout << "  and extending across the 123 sampling years.\n";
//
//  Compute the SVD.
//
  u = new double[m*m];
  s = new double[m*n];
  v = new double[n*n];

  r8mat_svd_linpack ( m, n, x, u, s, v );
//
//  Normalize the patterns so that each column has maximum entry 1.
//
  r8col_normalize_li ( n, n, v );
//
//  Reverse the row ordering.
//
  r8row_reverse ( n, n, v );
//
//  Plot the V modes.
//
  data.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    data << "  " << setw(4) << i
         << "  " << setw(14) << v[i+0*n]
         << "  " << setw(14) << v[i+1*n]
         << "  " << setw(14) << v[i+2*n]
         << "  " << setw(14) << v[i+3*n]
         << "  " << setw(14) << v[i+4*n]
         << "  " << setw(14) << v[i+5*n] << "\n";
  }
  data.close ( );
  cout << "  Created data file \"" << data_filename << "\".\n";

  command.open ( command_filename.c_str ( ) );
  command << "# " << command_filename << "\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < " << command_filename << "\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'vmode1.png'\n";
  command << "set xlabel 'Year'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Yearly Snowfall Mode 1'\n";
  command << "set grid\n";
  command << "set style data points\n";
  command << "plot '" << data_filename << "' using 1:2 with points lt 3 pt 3\n";

  command << "set output 'vmode2.png'\n";
  command << "set xlabel 'Year'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Yearly Snowfall Mode 2'\n";
  command << "set grid\n";
  command << "set style data points\n";
  command << "plot '" << data_filename << "' using 1:3 with points lt 3 pt 3\n";

  command << "set output 'vmode3.png'\n";
  command << "set xlabel 'Year'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Yearly Snowfall Mode 3'\n";
  command << "set grid\n";
  command << "set style data points\n";
  command << "plot '" << data_filename << "' using 1:4 with points lt 3 pt 3\n";

  command << "set output 'vmode4.png'\n";
  command << "set xlabel 'Year'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Yearly Snowfall Mode 4'\n";
  command << "set grid\n";
  command << "set style data points\n";
  command << "plot '" << data_filename << "' using 1:5 with points lt 3 pt 3\n";

  command << "set output 'vmode5.png'\n";
  command << "set xlabel 'Year'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Yearly Snowfall Mode 5'\n";
  command << "set grid\n";
  command << "set style data points\n";
  command << "plot '" << data_filename << "' using 1:6 with points lt 3 pt 3\n";

  command << "set output 'vmode6.png'\n";
  command << "set xlabel 'Year'\n";
  command << "set ylabel 'Snowfall'\n";
  command << "set title 'Yearly Snowfall Mode 6'\n";
  command << "set grid\n";
  command << "set style data points\n";
  command << "plot '" << data_filename << "' using 1:7 with points lt 3 pt 3\n";

  command << "quit\n";
  command.close ( );
  cout << "  Created command file '" << command_filename << "'.\n";

  delete [] s;
  delete [] u;
  delete [] v;

  return;
}
