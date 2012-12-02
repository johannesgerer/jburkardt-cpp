# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "poisson_simulation.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_SIMULATION_TEST tests POISSON_SIMULATION.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "POISSON_SIMULATION_TEST\n";
  cout << "  C++ version.\n";
  cout << "  Test the POISSON_SIMULATION library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "POISSON_SIMULATION_TEST\n";
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
//    TEST01 simulates waiting for a given number of events.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  int bin_num = 30;
  string command_filename;
  ofstream command;
  string data_filename;
  ofstream data;
  int event_num = 1000;
  int *f_bin;
  int i;
  int j;
  double lambda;
  int seed;
  double *t;
  double *w;
  double w_ave;
  double *w_bin;
  double w_max;
  double w_min;
  double width;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  POISSON_FIXED_EVENTS simulates a Poisson process\n";
  cout << "  until a given number of events have occurred.\n";
  cout << "\n";
  cout << "  Simulate a Poisson process, for which, on average,\n";
  cout << "  LAMBDA events occur per unit time.\n";
  cout << "  Run until you have observed EVENT_NUM events.\n";
 
  lambda = 0.5;
  seed = 123456789;

  cout << "\n";
  cout << "  LAMBDA = " << lambda << "\n";
  cout << "  EVENT_NUM = " << event_num << "\n";

  t = new double[event_num+1];
  w = new double[event_num+1];
  poisson_fixed_events ( lambda, event_num, seed, t, w );

  w_min = r8vec_min ( event_num + 1, w );
  w_max = r8vec_max ( event_num + 1, w );
  w_ave = r8vec_mean ( event_num + 1, w );

  cout << "\n";
  cout << "  Minimum wait = " << w_min << "\n";
  cout << "  Average wait = " << w_ave << "\n";
  cout << "  Maximum wait = " << w_max << "\n";

  cout << "\n";
  cout << " Count            Time            Wait\n";
  cout << "\n";
  for ( i = 0; i <= 5; i++ )
  {
    cout << "  " << i
         << "  " << t[i]
         << "  " << w[i] << "\n";
  }
  cout << "  ....  ..............  ..............\n";
  for ( i = event_num - 5; i <= event_num; i++ )
  {
    cout << "  " << i
         << "  " << t[i]
         << "  " << w[i] << "\n";
  }
//
//  Create the data file.
//
  data_filename = "poisson_timeline_data.txt";

  data.open ( data_filename.c_str ( ) );

  for ( i = 0; i <= event_num; i++ )
  {
    data << "  " << t[i]
         << "  " << i << "\n";
  }
  data.close ( );

  cout << " \n";
  cout << "  Data stored in \"" << data_filename << "\".\n";
//
//  Create the command file.
//
  command_filename = "poisson_timeline_commands.txt";

  command.open ( command_filename.c_str ( ) );

  command << "# poisson_timeline_commands.txt\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < poisson_timeline_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'poisson_timeline.png'\n";
  command << "set style data lines\n";
  command << "set xlabel 'Time'\n";
  command << "set ylabel 'Number of Events'\n";
  command << "set title 'Observation of Fixed Number of Poisson Events'\n";
  command << "set grid\n";
  command << "plot 'poisson_timeline_data.txt' using 1:2 lw 2\n";
  command << "quit\n";

  command.close ( );

  cout << "  Plot commands stored in \"" << command_filename << "\".\n";
//
//  Determine bin information.
//
  w_min = r8vec_min ( event_num + 1, w );
  w_max = r8vec_max ( event_num + 1, w );

  w_bin = r8vec_midspace_new ( bin_num, w_min, w_max );
  f_bin = new int[bin_num];

  for ( i = 0; i < bin_num; i++ )
  {
    f_bin[i] = 0;
  }

  for ( i = 0; i <= event_num; i++ )
  {
    j = 1 + ( int ) ( ( double ) ( bin_num ) * ( w[i] - w_min ) / ( w_max - w_min ) );
    j = i4_min ( j, bin_num );
    f_bin[j] = f_bin[j] + 1;
  }
//
//  Create the data file.
//
  data_filename = "poisson_times_data.txt";

  data.open ( data_filename.c_str ( ) );

  for ( i = 0; i < bin_num; i++ )
  {
    data << "  " << w_bin[i]
         << "  " << f_bin[i] << "\n";
  }
  data.close ( );

  cout << " \n";
  cout << "  Data stored in \"" << data_filename << "\".\n";
//
//  Create the command file.
//
  command_filename = "poisson_times_commands.txt";

  command.open ( command_filename.c_str ( ) );

  command << "# poisson_times_commands.txt\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < poisson_times_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'poisson_times.png'\n";
  command << "set xlabel 'Waiting Time'\n";
  command << "set ylabel 'Frequency'\n";
  command << "set title 'Waiting Times Observed Over Fixed Time'\n";
  command << "set grid\n";
  command << "set style fill solid\n";
  width = 0.85 * ( w_max - w_min ) / ( double ) ( bin_num );
  command << "plot 'poisson_times_data.txt' using 1:2:(" << width << ") with boxes\n";
  command << "quit\n";

  command.close ( );

  cout << "  Plot commands stored in \"" << command_filename << "\".\n";

  delete [] f_bin;
  delete [] t;
  delete [] w;
  delete [] w_bin;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 simulates waiting for a given length of time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  int bin_num = 30;
  string command_filename;
  ofstream command;
  string data_filename;
  ofstream data;
  int *f_bin;
  int i;
  double lambda;
  int *n;
  double *n_bin;
  double n_max;
  double n_mean;
  double n_min;
  double n_var;
  int seed;
  double t;
  int test;
  int test_num = 20000;
  double w;

  lambda = 0.5;
  t = 1000.0;
  seed = 123456789;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  POISSON_FIXED_EVENTS simulates a Poisson process\n";
  cout << "  counting the number of events that occur during\n";
  cout << "  a given time.\n";
  cout << "\n";
  cout << "  Simulate a Poisson process, for which, on average,\n";
  cout << "  LAMBDA events occur per unit time.\n";
  cout << "  Run for a total of " << t << " time units.\n";
  cout << "  LAMBDA = " << lambda << "\n";

  n = new int[test_num];

  for ( test = 0; test < test_num; test++ )
  {
    n[test] = poisson_fixed_time ( lambda, t, seed );
  }

  n_mean = i4vec_mean ( test_num, n );
  n_var = i4vec_variance ( test_num, n );
  cout << "\n";
  cout << "  Mean number of events = " << n_mean << "\n";
  cout << "  Variance = " << n_var << "\n";
  cout << "  STD = " << sqrt ( n_var ) << "\n";

  n_min = ( double ) ( i4vec_min ( test_num, n ) );
  n_max = ( double ) ( i4vec_max ( test_num, n ) );

  n_bin = r8vec_midspace_new ( bin_num, n_min, n_max );

  f_bin = new int[bin_num];
  for ( i = 0; i < bin_num; i++ )
  {
    f_bin[i] = 0;
  }
  for ( test = 0; test < test_num; test++ )
  {
    i = 1 + ( int ) ( ( double ) ( bin_num * ( n[test] - n_min ) ) 
      / ( double ) ( n_max - n_min ) );
    i = i4_min ( i, bin_num );
    f_bin[i] = f_bin[i] + 1;
  }
//
//  Create the data file.
//
  data_filename = "poisson_events_data.txt";

  data.open ( data_filename.c_str ( ) );

  for ( i = 0; i < bin_num; i++ )
  {
    data << "  " << n_bin[i]
         << "  " << f_bin[i] << "\n";
  }
  data.close ( );

  cout << " \n";
  cout << "  Data stored in \"" << data_filename << "\".\n";
//
//  Create the command file.
//
  command_filename = "poisson_events_commands.txt";

  command.open ( command_filename.c_str ( ) );

  command << "# poisson_events_commands.txt\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < poisson_events_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'poisson_events.png'\n";
  command << "set xlabel 'Number of Events'\n";
  command << "set ylabel 'Frequency'\n";
  command << "set title 'Number of Poisson Events Over Fixed Time'\n";
  command << "set grid\n";
  command << "set style fill solid\n";
  w = 0.85 * ( n_max - n_min ) / ( double ) ( bin_num );
  command << "plot 'poisson_events_data.txt' using 1:2:(" << w << ") with boxes\n";
  command << "quit\n";

  command.close ( );

  cout << "  Plot commands stored in \"" << command_filename << "\".\n";

  delete [] f_bin;
  delete [] n;
  delete [] n_bin;

  return;
}
