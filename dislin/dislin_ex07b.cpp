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
# define N 18

  char cbuf[80];
  int i;
  int icray[N] = { 30, 30, 30, 30, 30, 30, 100, 100, 100, 100,
                 100, 100, 170, 170, 170, 170, 170, 170 };
  float xray[N]  = { 1.0, 3.0, 8.0, 1.5, 9.0, 6.3, 5.8, 2.3, 8.1, 3.5,
                     2.2, 8.7, 9.2, 4.8, 3.4, 6.9, 7.5, 3.8 };
  float xwray[N];
  float yray[N]  = {5.0, 8.0, 3.5, 2.0, 7.0, 1.0, 4.3, 7.2, 6.0, 8.5,
                    4.1, 5.0, 7.3, 2.8, 1.6, 8.9, 9.5, 3.2 };
  float ywray[N];
  float z1ray[N] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  float z2ray[N] = { 4.0, 5.0, 3.0, 2.0, 3.5, 4.5, 2.0, 1.6, 3.8, 4.7,
                    2.1, 3.5, 1.9, 4.2, 4.9, 2.8, 3.6, 4.3 };

  cout << "\n";
  cout << "DISLIN_EX07B:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate of 3D bar graphs.\n";

  for ( i = 0; i < N; i++ )
  {
    xwray[i] = 0.5;
    ywray[i] = 0.5;
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
  setfil ( "dislin_ex07b.png" );
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
  axspos ( 200, 2600 );
  axslen ( 1800, 1800 );

  name ( "X-axis", "x" );
  name ( "Y-axis", "y" );
  name ( "Z-axis", "z" );

  titlin ( "3-D Bars / BARS3D", 3 );
  labl3d ( "hori" );

  graf3d ( 0.0, 10.0, 0.0, 2.0, 0.0, 10.0, 0.0, 2.0, 0.0, 5.0, 0.0, 1.0 );
  grid3d ( 1, 1, "bottom" );
  bars3d ( xray, yray, z1ray, z2ray, xwray, ywray, icray, N );

  legini ( cbuf, 3, 20 );
  legtit ( " ");
  legpos ( 1350, 1150 );
  leglin ( cbuf, "First", 1 );
  leglin ( cbuf, "Second", 2 );
  leglin ( cbuf, "Third", 3 );
  legend ( cbuf, 3);

  height ( 50 );
  title ( );
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX07B:\n";
  cout << "  Normal end of execution.\n";

  return 0;

# undef N
} 
