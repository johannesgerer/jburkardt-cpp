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
//    MAIN demonstrates the creation of a contour plot.
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
# define N 50

  float fpi = 3.14159/180.0;
  int i;
  int j;
  float step;
  float x;
  float xray[N];
  float y;
  float yray[N];
  float zlev;
  float zmat[N][N];

  cout << "\n";
  cout << "DISLIN_EX11:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the creation of a contour plot.\n";

  step = 360.0 / ( float ) ( N - 1 );

  for ( i = 0; i < N; i++ )
  {
    xray[i] = ( float ) i * step;
    yray[i] = ( float ) i * step;
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      x = xray[i] * fpi;
      y = yray[j] * fpi;    
      zmat[i][j] = 2.0 * sin(x) * sin(y);
    }
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
  setfil ( "dislin_ex11.png" );
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

  titlin ( "Contour Plot", 1 );
  titlin ( "F(X,Y) = 2 * SIN(X) * SIN(Y)", 3 );

  name ( "X-axis", "x" );
  name ( "Y-axis", "y" );

  intax ( );
  axspos ( 450, 2670 );
  graf ( 0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0 );

  height ( 30 );

  for ( i = 0; i < 9; i++ )
  {
    zlev = -2.+i*0.5;
    setclr ( ( i + 1 ) * 25 );
    if ( i == 4 )
    {
      labels ( "none", "contur" ); 
    }
    else
    {
      labels ( "float", "contur" );
    }

    contur ( xray, N, yray, N, (float *) zmat, zlev );
  }

  height ( 50 );
  color ( "fore" );
  title ( );
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX11:\n";
  cout << "  Normal end of execution.\n";

  return 0;

# undef N
}
