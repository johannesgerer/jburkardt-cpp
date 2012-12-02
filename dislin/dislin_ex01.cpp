# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "dislin.hpp"

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN demonstrates the use of the CURVE routine.
//
//  Modified:
//
//    09 April 2011
//
//  Reference:
//
//    Helmut Michels,
//    The Data Plotting Software DISLIN - version 10.4,
//    Shaker Media GmbH, January 2010,
//    ISBN13: 978-3-86858-517-9.
//
{
# define N 100

  double fpi = 3.1415926 / 180.0; 
  int i;
  double step;
  double x;
  float xray[N];
  float y1ray[N];
  float y2ray[N];

  cout << "\n";
  cout << "DISLIN_EX01:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the use of the CURVE routine, for\n";
  cout << "  plotting (X,Y) data.\n";

  step = 360.0 / ( double ) ( N - 1 );

  for ( i = 0; i < N; i++ )
  {
    xray[i] = i * step;
    x = xray[i] * fpi;
    y1ray[i] = sin ( x );
    y2ray[i] = cos ( x );
  }
//
//  Specify the format of the output file.
//
  metafl ( "png" );
//
//  Indicate that new data overwrites old data.
//
  filmod ( "delete" );
//
//  Specify the name of the output graphics file.
//
  setfil ( "dislin_ex01.png" );
//
//  Choose the page size and orientation.
//
  setpag ( "usal" );
//
//  For PNG output, reverse the default black background to white.
//
  scrmod ( "reverse" );
//
//  Open DISLIN.
//
  disini ( );
//
//  Plot a border around the page.
//
  pagera ( );
//
//  Use the COMPLEX font.
//
  complx ( );
  axspos ( 450, 1800 );
  axslen ( 2200, 1200 );

  name ( "X-axis", "x" );
  name ( "Y-axis", "y" );

  labdig ( -1, "x" );
  ticks ( 10, "xy" );

  titlin ( "Demonstration of CURVE", 1 );
  titlin ( "SIN(X), COS(X)", 3 );

  graf ( 0.0, 360.0, 0.0, 90.0, -1.0, 1.0, -1.0, 0.5 );
  title ( );

  color ( "red" );
  curve ( xray, y1ray, N );

  color ( "green" );
  curve ( xray, y2ray, N );

  color ( "fore" );
  dash ( );
  xaxgit ( );
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX01:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
