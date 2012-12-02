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
//    MAIN demonstrates the creation of 3D bar and pie graphs.
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
  char cbuf[80];
  int ic1ray[5]  = { 50, 150, 100, 200, 175 };
  int ic2ray[5]  = { 50, 150, 100, 200, 175 };
  float xray[5]  = { 2.0, 4.0, 6.0, 8.0, 10.0 };
  float y1ray[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
  float y2ray[5] = { 3.2, 1.5, 2.0, 1.0, 3.0 };

  cout << "\n";
  cout << "DISLIN_EX07:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the creation of 3D bar and pie graphs.\n";
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
  setfil ( "dislin_ex07.png" );
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

  titlin ( "3-D Bar Graph / 3-D Pie Chart", 2 );
  htitle ( 40 );

  shdpat ( 16 );
  axslen ( 1500, 1000 );
  axspos ( 300, 1400 );

  barwth ( 0.5);
  bartyp ( "3dvert" );
  labels ( "second", "bars" );
  labpos ( "outside", "bars" );
  labclr ( 255, "bars" );
  graf ( 0.0, 12.0, 0.0, 2.0, 0.0, 5.0, 0.0, 1.0 );
  title ( );
  color ( "red" );
  bars ( xray, y1ray, y2ray, 5 );
  endgrf ( );

  shdpat ( 16 );
  labels ( "data", "pie" );
  labclr ( 255, "pie" );
  chnpie ( "none" );
  pieclr ( ic1ray, ic2ray, 5 );
  pietyp ( "3d" );
  axspos ( 300, 2700 );
  piegrf ( cbuf, 0, y2ray , 5 );
//
//  Close DISLIN.
//      
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX07:\n";
  cout << "  Normal end of execution.\n";

  return 0;

}
