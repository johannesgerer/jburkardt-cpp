# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "curve_plot.h"

//****************************************************************************80

void curve_plot ( int n, double x[], double y[], string name )

//****************************************************************************80
//
//  Purpose:
//
//    CURVE_PLOT creates files that allow GNUPLOT to plot a curve.
//
//  Discussion:
//
//     Given N sets of X and Y data, and a name "NAME", this function
//     creates two files:
//       NAME_data.txt, containing the plot data;
//       NAME_commands.txt, containing commands to gnuplot.
//     Once these files are created, the command
//       gnuplot < NAME_commands.txt
//     should start up gnuplot, and cause it to create the file
//       NAME.png
//     containing a plot of the data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data points.
//
//    Input, double X[N], the X values.
//
//    Input, double Y[N], the Y values.
//
//    Input, string NAME, a name to use for the files to be created.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  string plot_filename;
//
//  Write the data file.
//
  data_filename = name + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    data_unit << x[i] << "  "
              << y[i] << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Plot data written to the file \"" << data_filename << "\".\n";
//
//  Write the command file.
//
  command_filename = name + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "set term png\n";
  plot_filename = name + ".png";
  command_unit << "set output \"" << plot_filename << "\"\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "unset key\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Y--->'\n";
  command_unit << "set timestamp\n";
  command_unit << "plot \"" << data_filename << "\" using 1:2 with lines lw 3\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Command data written to \"" << command_filename << "\".\n";

  return;
}
