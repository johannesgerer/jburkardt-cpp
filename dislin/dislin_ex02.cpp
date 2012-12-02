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
//    MAIN demonstrates the use of the POLAR routine.
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
# define M 10
# define N 300

  float a;
  float f = 3.1415926 / 180.0;
  int i;
  double step;
  float x2[M];
  float xray[N];
  float y2[M];
  float yray[N];
  
  cout << "\n";
  cout << "DISLIN_EX02:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the use of the POLAR routine, for\n";
  cout << "  plotting (R,Theta) data.\n";

  step = 360.0 / ( float ) ( N - 1 );

  for ( i = 0; i < N; i++ )
  { 
    a = ( float ) i * step * f;
    yray[i] = a;
    xray[i] = sin ( 5.0 * a );
  }

  for ( i = 0; i < M; i++ )
  {
    x2[i] = ( float ) ( i + 1 );
    y2[i] = ( float ) ( i + 1 );
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
  setfil ( "dislin_ex02.png" );
//
//  Choose the page size and orientation.
//
  setpag ( "usap" );
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
//  Use the HARDWARE font.
//
  hwfont ( );
  axspos ( 450, 1800 );

  titlin ( "Polar Plots", 2 );
  ticks  ( 3, "Y" );
  axends ( "NOENDS", "X" );
  labdig ( -1, "Y" );
  axslen ( 1000, 1000 );
  axsorg ( 1050, 900 );

  polar  ( 1.0, 0.0, 0.2, 0.0, 30.0 );
  curve  ( xray, yray, N );
  htitle ( 50 );
  title  ( );
  endgrf ( );

  labdig ( -1, "X" );
  axsorg ( 1050, 2250 );
  labtyp ( "VERT", "Y" );
  polar  ( 10.0, 0.0, 2.0, 0.0, 30.0 );
  barwth ( -5.0 );
  polcrv ( "FBARS" );
  curve  ( x2, y2, M );
//
//  Close DISLIN.
//       
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX02:\n";
  cout << "  Normal end of execution.\n";

  return 0;
# undef M
# undef N
}
