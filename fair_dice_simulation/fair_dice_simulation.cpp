# include <cstdlib>
# include <iostream>
# include <fstream>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
int random_int ( int a, int b );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    FAIR_DICE_SIMULATION simulates N throws of two fair dice.
//
//  Usage:
//
//    fair_dice n
//
//    where 
//
//    * n is the number of times the dice should be thrown.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Command line, int N, the number of times the dice are thrown.
//
{
  ofstream command;
  string command_filename = "fair_dice_commands.txt";
  ofstream data;
  string data_filename = "fair_dice_data.txt";
  int die1;
  int die2;
  int i;
  int n;
  int seed;
  int score;
  int score_count[13];

  timestamp ( );
  cout << "\n";
  cout << "FAIR_DICE_SIMULATION:\n";
  cout << "  C++ version\n";
  cout << "  Simulate N throws of a pair of fair dice.\n";

  if ( 1 < argc )
  {
    n = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter N, the number of times the dice are thrown: ";
    cin >> n;
  }
//
//  Initialize the random number generator.
//
  seed = time ( 0 );
  srand ( seed );
//
//  For convenience, include slots for 0 and 1, even though they never occur.
//
  for ( i = 0; i <= 12; i++ )
  {
    score_count[i] = 0;
  }
//
//  Roll N times.
//
  for ( i = 1; i <= n; i++ )
  {
    die1 = random_int ( 1, 6 );
    die2 = random_int ( 1, 6 );
    score = die1 + die2;
    score_count[score] = score_count[score] + 1;
  }
//
//  Write the graphics data file.
//
  data.open ( data_filename.c_str ( ) );
  for ( score = 2; score <= 12; score++ )
  {
    data << "  " << score << "  " << score_count[score] << "\n";
  }
  data.close ( );
  cout << "\n";
  cout << "  Created the graphics data file \"" << data_filename << "\".\n";
//
//  Write the graphics command file.
//
  command.open ( command_filename.c_str ( ) );
  command << "# " << command_filename << "\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < " << command_filename << "\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'fair_dice.png'\n";
  command << "set xlabel 'Score'\n";
  command << "set ylabel 'Frequency'\n";
  command << "set title 'Score frequency for a pair of fair dice'\n";
  command << "set grid\n";
  command << "set style fill solid\n";
  command << "set yrange [0:*]\n";
  command << "set timestamp\n";
  command << "plot 'fair_dice_data.txt' using 1:2:(0.90):xtic(3) with boxes\n";
  command << "quit\n";
  command.close ( );
  cout << "  Created the graphics command file \"" << command_filename << "\".\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "FAIR_DICE_SIMULATION:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int random_int ( int a, int b )

//****************************************************************************80
//
//  Purpose:
//
//   RANDOM_INT returns a random integer between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, B, the range of integers.
//
//    Output, int RANDOM_INT, the random integer.
//
{
  int range;
  int value;
//
//  If we want integers between A and B, there are actually
//  B - A + 1 values.
//
  range = b - a + 1;

  value = a + rand ( ) % range;

  return value;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
