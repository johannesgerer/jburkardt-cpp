# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

# include "gnuplot_i.hpp"

using namespace std;

# define SLEEP_LGTH  1

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ANIM.
//
//  Modified:
//
//    24 June 2011
//
{
  gnuplot_ctrl *h1;
  double phase;
 
  cout << "\n";
  cout << "ANIM:\n";
  cout << "  C++ version\n";
  cout << "  Demonstrate how a running C++ program can create plots\n";
  cout << "  during execution by calling gnuplot, using the\n";
  cout << "  gnuplot_i interface program.\n";
//
//  Open a GNUPLOT process.
//
  h1 = gnuplot_init ( );

  if ( h1 == NULL )
  {
    cout << "\n";
    cout << "ANIM - Fatal error!\n";
    cout << "  The gnuplot command is not in your path.\n"; 
    exit ( 1 );
  }

  for ( phase = 0.1; phase < 10.0; phase = phase + 0.1 )
  {
    gnuplot_resetplot ( h1 );
    gnuplot_cmd ( h1, "plot sin(x+%g)", phase );
  }

  for ( phase = 10.0; 0.0 <= phase; phase = phase - 0.1 ) 
  {
    gnuplot_resetplot ( h1 );
    gnuplot_cmd ( h1, "plot sin(x+%g)", phase );
  }
//
//  Close the GNUPLOT process.
//   
  gnuplot_close ( h1 );
//
//  Terminate.
//
  cout << "\n";
  cout << "ANIM:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}

