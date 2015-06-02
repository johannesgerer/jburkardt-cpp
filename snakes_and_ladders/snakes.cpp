# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

# include "snakes.hpp"

//****************************************************************************80

double *snakes ( )

//****************************************************************************80
//
//  Purpose:
//
//    SNAKES sets up the Snakes and Ladders matrix.
//
//  Discussion:
//
//    Snakes and Ladders, also known as Chutes and Ladders, is a game
//    played on a 10x10 board of 100 squares.  A player can be said to
//    start at square 0, that is, off the board.  The player repeatedly
//    rolls a die, and advances between 1 and 6 steps accordingly.
//    The game is won when the player reaches square 100.  In some versions,
//    the player must reach 100 by exact die count, forfeiting the move
//    if 100 is exceeded; in others, reaching or exceeding 100 counts as
//    a win.
//
//    Play is complicated by the existence of snakes and ladders.  On
//    landing on a square that begins a snake or ladder, the player is
//    immediately tranported to another square, which will be lower for
//    a snake, or higher for a ladder.
//
//    Typically, several players play, alternating turns.
//
//    Given a vector V(0:100) which is initially all zero except for the
//    first entry, the matrix-vector product A'*V represents the probabilities
//    that a player starting on square 0 will be on any given square after one
//    roll.  Correspondingly, (A')^2*V considers two moves, and so on.  Thus,
//    repeatedly multiplying by A' reveals the probability distribution 
//    associated with the likelihood of occupying any particular square at a 
//    given turn in the game.  
//
//    There is a single eigenvalue of value 1, whose corresponding eigenvector
//    is all zero except for a final entry of 1, representing a player who
//    has reached square 100.  All other eigenvalues have norm less than 1,
//    corresponding to the fact that there are no other long term steady
//    states or cycles in the game.
//
//    Note that no two descriptions of the Snakes and Ladders board seem to
//    agree.  This is the board described by Nick Berry.  The board described 
//    by Higham and Higham is close to this one, but differs in the description 
//    of two of the jumps.
//
//    While most commentators elect to move immediately from a snake mouth or
//    ladder foot, I have decide there are reasons to treat the game in such a
//    way that when you land on a ladder foot or snake mouth, you stay there
//    as though you had landed on an ordinary square; the difference arises on
//    your next turn, when, instead of rolling a die, you move up the ladder
//    or down the snake.  This allows the player to "register" a stop at the
//    given square, may be suitable for certain applications, and makes for
//    a transition matrix whose structure is more obvious to understand.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Steve Althoen, Larry King, Kenneth Schilling,
//    How long is a game of Snakes and Ladders?,
//    The Mathematical Gazette,
//    Volume 77, Number 478, March 1993, pages 71-76.
//
//    Nick Berry,
//    Mathematical Analysis of Chutes and Ladders,
//    http://www.datagenetics.com/blog/november12011/index.html
//
//    Desmond Higham, Nicholas Higham,
//    MATLAB Guide,
//    SIAM, 2005,
//    ISBN13: 9780898717891.
//
//  Parameters:
//
//    Output, double SNAKES[101*101], the matrix.
//
{
  double *a;
  int d;
  int i;
  int j;
  int j1;
  int j2;
  int *jump;
  int k;

  a = new double[101*101];

  jump = new int[101];

  for ( i = 0; i < 101; i++ )
  {
    jump[i] = i;
  }

  jump[ 1] =  38;
  jump[ 4] =  14;
  jump[ 9] =  31;
  jump[16] =   6;
  jump[21] =  42;
  jump[28] =  84;
  jump[36] =  44;
  jump[48] =  26;
  jump[49] =  11;
  jump[51] =  67;
  jump[56] =  53;
  jump[62] =  19;
  jump[64] =  60;
  jump[71] =  91;
  jump[80] = 100;
  jump[87] =  24;
  jump[93] =  73;
  jump[95] =  75;
  jump[98] =  78;

  for ( j = 0; j < 101; j++ )
  {
    for ( i = 0; i < 101; i++ )
    {
      a[i+j*101] = 0.0;
    }
  }
//
//  A(I,J) represents the probablity that a dice roll will take you from
//  square I to square J.
//
//  Starting in square I...
//
  for ( i = 0; i < 101; i++ )
  {
//
//  If I is a snake or ladder, go to the next spot.
//
    if ( i != jump[i] )
    {
      j = jump[i];
      a[i+j*101] = 1.0;
    }
//
//  Otherwise, roll a die
//
    else
    {
      for ( d = 1; d <= 6; d++ )
      {
//
//  so theoretically, our new location J will be I + D,
//
        j = i + d;
//
//  but if J is greater than 100, move us back to J,
//
        if ( 100 < j )
        {
          j = 100;
        }
  
        a[i+j*101] = a[i+j*101] + 1.0 / 6.0;
      }
    }
  }

  delete [] jump;

  return a;
}
//****************************************************************************80

void spy_ge ( int m, int n, double a[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    SPY_GE plots a sparsity pattern for a general storage (GE) matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns
//    in the matrix.
//
//    Input, double A[M*N], the matrix.
//
//    Input, string HEADER, the name to be used for the
//    title of the plot, and as part of the names of the data, command
//    and plot files.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  int j;
  int nz_num;
  string png_filename;
//
//  Create data file.
//
  data_filename = header + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  nz_num = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] != 0.0 )
      {
        data_unit << j << "  "
                  << i << "\n";
        nz_num = nz_num + 1;
      }
    }
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created sparsity data file '" << data_filename << "'\n";
//
//  Create command file.
//
  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "unset key\n";
  command_unit << "set term png\n";

  png_filename = header + ".png";
  command_unit << "set output '" << png_filename << "'\n";
  command_unit << "set size ratio -1\n";
  command_unit << "set xlabel '<--- J --->'\n";
  command_unit << "set ylabel '<--- I --->'\n";
  command_unit << "set title '" << nz_num << " nonzeros for \"" 
               << header << "\"'\n";
  command_unit << "set timestamp\n";
  command_unit << "plot [y=0:" << n - 1 << "] [x="
               << m - 1 << ":0] '"
               << data_filename << "' with points pt 5\n";
 
  command_unit.close ( );
  cout << "  Created graphics command file '" << command_filename << "'\n";

  return;
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
