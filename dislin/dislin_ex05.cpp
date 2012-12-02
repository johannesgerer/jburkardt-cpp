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
//    MAIN demonstrates the use of bar graphs.
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
  static char cbuf[25];
  static char *ctit = "Bar Graphs (BARS)";
  int i;
  int nya = 2700;
  static float x[9]  = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
  static float y[9]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  static float y1[9] = { 1.0, 1.5, 2.5, 1.3, 2.0, 1.2, 0.7, 1.4, 1.1 };
  static float y2[9] = { 2.0, 2.7, 3.5, 2.1, 3.2, 1.9, 2.0, 2.3, 1.8 };
  static float y3[9] = { 4.0, 3.5, 4.5, 3.7, 4.0, 2.9, 3.0, 3.2, 2.6 };

  cout << "\n";
  cout << "DISLIN_EX05:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the display of data in bar graphs.\n";
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
  setfil ( "dislin_ex05.png" );
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
//  Use the COMPLEX font.
//
  complx ( );
  ticks ( 1, "x" );
  intax ( );
  axslen ( 1600, 700 );
  titlin ( ctit, 3 );

  legini ( cbuf, 3, 8 );
  leglin ( cbuf, "FIRST", 1 );
  leglin ( cbuf, "SECOND", 2 );
  leglin ( cbuf, "THIRD", 3 );
  legtit ( " " );
  shdpat ( 5L );

  for ( i = 1; i <= 3; i++ )
  {
     if ( 1 < i )
     {
       labels ( "none", "x" );
     }
    axspos ( 300, nya-(i-1)*800 );

    graf ( 0.0, 10.0, 0.0, 1.0, 0.0, 5.0, 0.0, 1.0 );

    if ( i == 1 )
    {
      bargrp ( 3, 0.15 );
      color ( "red" );
      bars ( x, y, y1, 9 );
      color ( "green" );
      bars ( x, y, y2, 9 );
      color ( "blue" );
      bars ( x, y, y3, 9 );
      color ( "fore" );
      reset ( "bargrp" );
    }
    else if ( i == 2 )
    {
      height ( 30 );
      labels ( "delta", "bars" );
      labpos ( "center", "bars" );
      color ( "red" );
      bars ( x, y, y1, 9 );
      color ( "green" );
      bars ( x, y1, y2, 9 );
      color ( "blue" );
      bars ( x, y2, y3, 9 );
      color ( "fore" );
      reset ( "height" ); 
    }
    else if ( i == 3 )
    {
      labels ( "second", "bars" );
      labpos ( "outside", "bars" );
      color ( "red" );
      bars ( x, y, y1, 9 );
      color ( "fore" );
    }

    if ( i < 3 )
    {
      legend ( cbuf, 7 );
    }
    else
    {
      height ( 50 );
      title ( );
    }

    endgrf ( );
  }
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX05:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
