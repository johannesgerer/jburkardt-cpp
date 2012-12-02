# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <string>

using namespace std;

int main ( int argc, char *argv[] );
void i4_to_van_der_corput_sequence ( int seed, int base, int n, double r[] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for VAN_DER_CORPUT_DATASET.
//
//  Discussion:
//
//    VAN_DER_CORPUT_DATASET generates a van der Corput dataset.
//
//    This program is meant to be used interactively.  It's also
//    possible to prepare a simple input file beforehand and use it
//    in batch mode.
//
//  Usage:
//
//    van_der_corput_dataset base seed n
//
//    where
//
//    * BASE, the base of the sequence;
//    * SEED, the index of the first element to be computed;
//    * N, the number of points to generate.
//
//    The program generates the data and writes it to the file
//
//      van_der_corput_BASE_SEED_N.txt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 December 2009
//
//  Author:
//
//    John Burkardt
//
{
  int base;
  char output_filename[255];
  int i;
  int n;
  double *r;
  int seed;

  timestamp ( );

  cout << "\n";
  cout << "VAN_DER_CORPUT_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Generate a van der Corput dataset.\n";
  cout << "\n";
  cout << "  The program requests input values from the user:\n";
  cout << "\n";
  cout << "  * BASE, the base,\n";
  cout << "  * SEED, a positive integer.\n";
  cout << "  * N, the number of points to generate,\n";
  cout << "\n";
  cout << "  The program generates the data, writes it to the file\n";
  cout << "\n";
  cout << "    van_der_corput_BASE_SEED_N.txt\n";
  cout << "\n";
//
//  Get the BASE.
//
  if ( 1 < argc )
  {
    base = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter BASE, the van der Corput base,\n";
    cout << "  which is often a prime number:\n";
    cout << "  (Try '2' if you don't have a preference.)\n";
    cin >> base;
  }

  cout << "\n";
  cout << "  BASE = " << base << "\n";
//
//  Get SEED.
//
  if ( 2 < argc )
  {
    seed = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter SEED, which is 0 or a positive integer.\n";
    cout << "  (Try '0' or '1' if you don't have a preference.)\n";
    cin >> seed;
  }

  cout << "  SEED = " << seed << "\n";
//
//  Get N.
//
  if ( 3 < argc )
  {
    n = atoi ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter N, the number of points to generate:\n";
    cout << "  (Try '25' if you don't have a preference.)\n";
    cin >> n;
  }

  cout << "  N = " << n << "\n";
//
//  Compute the data.
//
  r = new double[n];

  i4_to_van_der_corput_sequence ( seed, base, n, r );
//
//  Write the data to a file.
//
  sprintf ( output_filename, "van_der_corput_%d_%d_%d.txt", base, seed, n );

  r8mat_write ( output_filename, 1, n, r );

  cout << "\n";
  cout << "  The data was written to the file \"" 
    << output_filename << "\".\n";
//
//  Terminate.
//
  delete [] r;

  cout << "\n";
  cout << "VAN_DER_CORPUT_DATASET:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void i4_to_van_der_corput_sequence ( int seed, int base, int n, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_VAN_DER_CORPUT_SEQUENCE: next N elements of a van der Corput sequence.
//
//  Modified:
//
//    27 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J H Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, pages 84-90, 1960.
//
//    J G van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Input, int SEED, the index of the first desired element.
//    SEED should be nonnegative.
//    SEED = 0 is allowed, and returns R = 0.
//
//    Input, int BASE, the van der Corput base, which is typically a
//    prime number.
//
//    Input, int N, the number of elements desired.
//
//    Output, double R[N], the SEED-th through 
//    (SEED+N-1)-th elements of the van der Corput sequence for 
//    the given base.
//
{
  double base_inv;
  int digit;
  int i;
  int seed2;

  if ( base <= 1 ) 
  {
    cout << "\n";
    cout << "I4_TO_VAN_DER_CORPUT_SEQUENCE - Fatal error!\n";
    cout << "  The input base BASE is <= 1!\n";
    cout << "  BASE = " << base << "\n";
    exit ( 1 );
  }

  if ( seed < 0 ) 
  {
    cout << "\n";
    cout << "I4_TO_VAN_DER_CORPUT_SEQUENCE - Fatal error!\n";
    cout << "  The input base SEED is < 0!\n";
    cout << "  SEED = " << seed << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    r[i] = 0.0;
    seed2 = seed + i;
    base_inv = 1.0E+00 / ( ( double ) base );

    while ( seed2 != 0 )
    {
      digit = seed2 % base;
      r[i] = r[i] + ( ( double ) digit ) * base_inv;
      base_inv = base_inv / ( ( double ) base );
      seed2 = seed2 / base;
    }
  }

  return;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
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
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    04 October 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

