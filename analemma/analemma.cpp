# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <ctime>
# include <fstream>

using namespace std;

int main ( int argc, char **argv );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char **argv )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ANALEMMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    Original C version by Brian Tung.
//    C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double ECC, the orbital eccentricity.
//
//    Local, double LON, the longitude of the perihelion in radians.
//
//    Local, double OBLIQ, the obliquity in radians.
//
{
  char c;
  string command_filename = "analemma_commands.txt";
  ofstream command_unit;
  string data_filename = "analemma_data.txt";
  ofstream data_unit;
  double days = 365.242;
  double dec;
  double degrees = ( 3.141592653589793 / 180.0 );
  double ecc = 0.01671;
  double eot;
  double f;
  double lon = 1.347;
  double obliq = 0.4091;
  extern char *optarg;
  double pi = 3.141592653589793;
  double t;
  double tau;
  double theta;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;
  double z1;
  double z2;
  double z3;

  timestamp ( );
  cout << "\n";
  cout << "ANALEMMA\n";
  cout << "  C++ version\n";
  cout << "  Compute and plot the analemma, equation of time, and declination.\n";
  cout << "  This program is based on a C program by Brian Tung.\n";
/* 
  Parse the arguments 
*/
  while ( ( c = getopt ( argc, argv, "e:l:o:h" ) ) >= 0)
  {
    switch ( c ) 
    {
      case 'e':
        ecc = atof ( optarg );
        break;
      case 'l':
        lon = atof ( optarg ) * degrees;
        break;
      case 'o':
        obliq = atof ( optarg ) * degrees;
        break;
      default:
        cerr << "Usage: analemma [options]\n";
        cerr << "    -e <ecc>    eccentricity ";
        cerr << "(default value: " << ecc << ")\n";
        cerr << "    -l <lon>    longitude of perihelion in deg ";
        cerr << "(default value: " << lon / degrees << "\n";
        cerr << "    -o <obliq>  axial obliquity in deg ";
        cerr << "(default value: " << obliq / degrees << "\n";
        cerr << "    -h          print this page\n";
        exit ( 0 );
    }
  }
//
//  Compute the data.
//
  data_unit.open ( data_filename.c_str ( ) );

  for ( f = 0.0; f <= 1.0; f = f + 0.0001 ) 
  {
    tau = 2.0 * pi * f;
//
//  Set theta to the current longitude. 
//
    theta = atan2 ( sqrt ( 1.0 - ecc * ecc ) * sin ( tau ), cos ( tau ) - ecc );
//
//  Rotate clockwise in XY plane by theta, corrected by lon.
//
    x1 = cos ( theta - ( lon - pi / 2.0 ) );
    y1 = sin ( theta - ( lon - pi / 2.0 ) );
    z1 = 0.0;
//
//  Rotate counter-clockwise in XZ plane by obliq.
//
    x2 = cos ( obliq ) * x1 + sin ( obliq ) * z1;
    y2 = y1;
    z2 = - sin ( obliq ) * x1 + cos ( obliq ) * z1;
// 
//  Set t equal to real time from tau and
//  rotate counter-clockwise by t, corrected by lon 
//
    t = tau - ecc * sin ( tau );
    x3 =   cos ( t - ( lon - pi / 2.0 ) ) * x2 + sin ( t - ( lon - pi / 2.0 ) ) * y2;
    y3 = - sin ( t - ( lon - pi / 2.0 ) ) * x2 + cos ( t - ( lon - pi / 2.0 ) ) * y2;
    z3 = z2;

    eot = - atan2 ( y3, x3 ) * 4.0 / degrees * days / ( days + 1.0 );
    dec = asin ( z3 ) / degrees;
// 
//  Print results in minutes early/late and degrees north/south 
//
    data_unit << "  " << t / ( 2.0 * pi )
              << "  " << eot
              << "  " << dec << "\n";
  }

  data_unit.close ( );

  cout << "\n";
  cout << "  Created data file \"" << data_filename << "\".\n";
//
//  Create the command file.
//
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "set term png\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "unset key\n";
  command_unit << "set output \"eot.png\"\n";
  command_unit << "set xlabel '<---Normalized Date--->'\n";
  command_unit << "set ylabel '<---Minutes Early/Late--->'\n";
  command_unit << "set title 'The equation of time'\n";
  command_unit << "plot '" << data_filename << "' using 1:2 with lines\n";
  command_unit << "set output \"declination.png\"\n";
  command_unit << "set xlabel '<---Normalized Date--->'\n";
  command_unit << "set ylabel '<---Degrees North/South--->'\n";
  command_unit << "set title 'Declination'\n";
  command_unit << "plot '" << data_filename << "' using 1:3 with lines\n";
  command_unit << "set output \"analemma.png\"\n";
  command_unit << "set xlabel '<---Minutes Early/Late--->'\n";
  command_unit << "set ylabel '<---Degrees North/South--->'\n";
  command_unit << "set title 'The analemma'\n";
  command_unit << "plot '" << data_filename << "' using 2:3 with lines\n";
  command_unit << "quit\n";
  command_unit << "\n";

  command_unit.close ( );

  cout << "  Created command file \"" << command_filename << "\".\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "ANALEMMA\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
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
