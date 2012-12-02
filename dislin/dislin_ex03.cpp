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
//    MAIN demonstrates the use of the SYMBOL routine.
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
  static char cstr[3];
  static char ctit[] = "Symbols";
  int i;
  int nl;
  int nxp;
  int ny;

  cout << "\n";
  cout << "DISLIN_EX03:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the use of the SYMBOL routine, for\n";
  cout << "  the generation of special symbols and fonts.\n";
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
  setfil ( "dislin_ex03.png" );
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

  height ( 60 );
  nl = nlmess ( ctit );
  messag ( ctit, ( 2100 - nl ) / 2, 200 );

  height ( 50 );
  hsymbl ( 120 );

  ny = 150;

  for ( i = 0; i < 24; i++ )
  {
    if ( ( i % 4 ) == 0 ) 
    {
      ny  = ny + 400;
      nxp = 550;
    }
    else
    {
      nxp = nxp + 350;
    }

    sprintf ( cstr, "%d", i ); 
    nl = nlmess ( cstr ) / 2;
    messag ( cstr, nxp-nl ,ny+150 );
    symbol ( i, nxp, ny );
  }
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX03:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
