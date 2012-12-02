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
//    MAIN demonstrates the use of shade patterns.
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
  static char *ctit = "Shading Patterns (AREAF)";
  int i;
  int iclr;
  int ii;
  static int ix[4] = { 0, 300, 300, 0 };
  int ixp[4];
  static int iy[4] = {0, 0, 400, 400 };
  int iyp[4];
  int j;
  int k;
  int nl;
  int nx;
  int nx0 = 335;
  int ny;
  int ny0 = 350;

  cout << "\n";
  cout << "DISLIN_EX08:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the use of shade patterns.\n";
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
  setfil ( "dislin_ex08.png" );
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
  setvlt ( "small" );

  height ( 50 );
  nl = nlmess ( ctit );
  nx = ( 2970 - nl ) / 2 ;
  messag ( ctit, nx, 200 );

  iclr = 0;
  for ( i = 0; i < 3; i++ )
  {
    ny = ny0 + i * 600;
    for ( j = 0; j < 6; j++ )
    {
      nx = nx0 + j * 400;
      ii = i * 6 + j;
      shdpat ( (long) ii );
      sprintf ( cstr, "%d", ii );

      iclr = iclr % 16;
      iclr = iclr + 1;
      setclr ( iclr );

      for ( k = 0; k < 4; k++ )
      {
        ixp[k] = ix[k] + nx;
        iyp[k] = iy[k] + ny;
      }
      areaf ( ixp, iyp, 4 );

      nl = nlmess ( cstr );
      nx = nx + ( 300 - nl ) / 2;
      messag ( cstr, nx, ny+460 );
    }
  }
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX08:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
