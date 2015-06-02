# include <cstdlib>
# include <iostream>
# include <fstream>
# include <sstream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "vandermonde_interp_1d.hpp"
# include "condition.hpp"
# include "qr_solve.hpp"
# include "test_interp.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob );
void test02 ( int prob );
string i4_to_string ( int i4 );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for VANDERMONDE_INTERP_1D_PRB.
//
//  Discussion:
//
//    VANDERMONDE_INTERP_1D_PRB tests the VANDERMONDE_INTERP_1D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "VANDERMONDE_INTERP_1D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the VANDERMONDE_INTERP_1D library.\n";
  cout << "  The QR_SOLVE library is needed.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  This test needs the CONDITION library.\n";
  cout << "  This test needs the TEST_INTERP library.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    test01 ( prob );
  }

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    test02 ( prob );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "VANDERMONDE_INTERP_1D_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests VANDERMONDE_INTERP_1D_MATRIX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int debug = 0;
  double *c;
  double condition;
  int i;
  double int_error;
  double ld;
  double li;
  int m;
  int nd;
  int ni;
  double *xd;
  double *xi;
  double xmax;
  double xmin;
  double *xy;
  double *yd;
  double *yi;
  double ymax;
  double ymin;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Interpolate data from TEST_INTERP problem #" << prob << "\n";

  nd = p00_data_num ( prob );
  cout << "  Number of data points = " << nd << "\n";

  xy = p00_data ( prob, 2, nd );
  
  if ( debug )
  {
    r8mat_transpose_print ( 2, nd, xy, "  Data array:" );
  }

  xd = new double[nd];
  yd = new double[nd];

  for ( i = 0; i < nd; i++ )
  {
    xd[i] = xy[0+i*2];
    yd[i] = xy[1+i*2];
  }
//
//  Choose the degree of the polynomial to be ND - 1.
//
  m = nd - 1;
//
//  Compute Vandermonde matrix and get condition number.
//
  a = vandermonde_interp_1d_matrix ( nd, xd );

  condition = condition_hager ( nd, a );

  cout << "\n";
  cout << "  Condition of Vandermonde matrix is " << condition << "\n";
//
//  Solve linear system.
//
  c = qr_solve ( nd, nd, a, yd );
//
//  #1:  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = r8poly_value ( m, c, ni, xi );

  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xi;
  delete [] yi;
//
//  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
//  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
//  (YMAX-YMIN).
//
  xmin = r8vec_min ( nd, xd );
  xmax = r8vec_max ( nd, xd );
  ymin = r8vec_min ( nd, yd );
  ymax = r8vec_max ( nd, yd );

  ni = 501;
  xi = r8vec_linspace_new ( ni, xmin, xmax );
  yi = r8poly_value ( m, c, ni, xi );

  ld = 0.0;
  for ( i = 0; i < nd - 1; i++ )
  {
    ld = ld + sqrt ( pow ( ( xd[i+1] - xd[i] ) / ( xmax - xmin ), 2 )
                   + pow ( ( yd[i+1] - yd[i] ) / ( ymax - ymin ), 2 ) ); 
  }

  li = 0.0;
  for ( i = 0; i < ni - 1; i++ )
  {
    li = li + sqrt ( pow ( ( xi[i+1] - xi[i] ) / ( xmax - xmin ), 2 )
                   + pow ( ( yi[i+1] - yi[i] ) / ( ymax - ymin ), 2 ) );
  }

  cout << "\n";
  cout << "  Normalized length of piecewise linear interpolant = " << ld << "\n";
  cout << "  Normalized length of polynomial interpolant       = " << li << "\n";

  delete [] a;
  delete [] c;
  delete [] xd;
  delete [] xi;
  delete [] xy;
  delete [] yd;
  delete [] yi;

  return;
}
//****************************************************************************80

void test02 ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests VANDERMONDE_INTERP_1D_MATRIX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
{
  double *a;
  double *c;
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  string interp_filename;
  ofstream interp_unit;
  int j;
  int nd;
  int ni;
  string output_filename;
  string title;
  double *xd;
  double *xi;
  double xmax;
  double xmin;
  double *xy;
  double *yd;
  double *yi;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  VANDERMONDE_INTERP_1D_MATRIX sets the Vandermonde linear system\n";
  cout << "  for the interpolating polynomial.\n";
  cout << "  Interpolate data from TEST_INTERP problem #" << prob << "\n";

  nd = p00_data_num ( prob );
  cout << "  Number of data points = " << nd << "\n";

  xy = p00_data ( prob, 2, nd );
  
  r8mat_transpose_print ( 2, nd, xy, "  Data array:" );

  xd = new double[nd];
  yd = new double[nd];

  for ( i = 0; i < nd; i++ )
  {
    xd[i] = xy[0+2*i];
    yd[i] = xy[1+2*i];
  }
//
//  Compute Vandermonde matrix and get condition number.
//
  a = vandermonde_interp_1d_matrix ( nd, xd );
//
//  Solve linear system.
//
  c = qr_solve ( nd, nd, a, yd );
//
//  Create data file.
//
  data_filename = "data" + i4_to_string ( prob ) + ".txt";
  data_unit.open ( data_filename.c_str ( ) );
  for ( j = 0; j < nd; j++ )
  {
    data_unit << "  " << xd[j]
              << "  " << yd[j] << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created graphics data file \"" << data_filename << "\".\n";
//
//  Create interp file.
//
  ni = 501;
  xmin = r8vec_min ( nd, xd );
  xmax = r8vec_max ( nd, xd );
  xi = r8vec_linspace_new ( ni, xmin, xmax );
  yi = r8poly_value ( nd - 1, c, ni, xi );

  interp_filename = "interp" + i4_to_string ( prob ) + ".txt";
  interp_unit.open ( interp_filename.c_str ( ) );
  for ( j = 0; j < ni; j++ )
  {
    interp_unit << "  " << xi[j]
                << "  " << yi[j] << "\n";
  }
  interp_unit.close ( );
  cout << "  Created graphics interp file \"" << interp_filename << "\".\n";
//
//  Plot the data and the interpolant.
//
  command_filename = "commands" + i4_to_string ( prob ) + ".txt";
  command_unit.open ( command_filename.c_str ( ) );

  output_filename = "plot" + i4_to_string ( prob ) + ".png";

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output '" << output_filename << "'\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Y--->'\n";
  command_unit << "set title 'Data versus Vandermonde polynomial interpolant'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2 with points pt 7 ps 2 lc rgb 'blue',\\\n";
  command_unit << "     '" << interp_filename 
               << "' using 1:2 lw 3 linecolor rgb 'red'\n";

  command_unit.close ( );
  cout << "  Created graphics command file \"" << command_filename << "\".\n";
//
//  Free memory.
//
  delete [] xd;
  delete [] xi;
  delete [] xy;
  delete [] yd;
  delete [] yi;

  return;
}
//****************************************************************************80

string i4_to_string ( int i4 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  ostringstream fred;
  string value;

  fred << i4;

  value = fred.str ( );

  return value;
}
