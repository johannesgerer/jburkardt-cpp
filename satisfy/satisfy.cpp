# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
int circuit_value ( int n, int bvec[] );
void i4_to_bvec ( int i4, int n, int bvec[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SATISFY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Michael Quinn,
//    Parallel Programming in C with MPI and OpenMP,
//    McGraw-Hill, 2004,
//    ISBN13: 978-0071232654,
//    LC: QA76.73.C15.Q55.
//
{
# define N 23

  int bvec[N];
  int i;
  int ihi;
  int j;
  int n = N;
  int solution_num;
  int value;

  cout << "\n";
  timestamp ( );
  cout << "\n";
  cout << "SATISFY\n";
  cout << "  C++ version\n" ;
  cout << "  We have a logical function of N logical arguments.\n";
  cout << "  We do an exhaustive search of all 2^N possibilities,\n";
  cout << "  seeking those inputs that make the function TRUE.\n";
//
//  Compute the number of binary vectors to check.
//
  ihi = 1;
  for ( i = 1; i <= n; i++ )
  {
    ihi = ihi * 2;
  }
  cout << "\n";
  cout << "  The number of logical variables is N = " << n << "\n";;
  cout << "  The number of input vectors to check is " << ihi << "\n";
  cout << "\n";
  cout << "   #       Index    ---------Input Values------------------------\n";
  cout << "\n";
//
//  Check every possible input vector.
//
  solution_num = 0;

  for ( i = 0; i < ihi; i++ )
  {
    i4_to_bvec ( i, n, bvec );

    value = circuit_value ( n, bvec );

    if ( value == 1 )
    {
      solution_num = solution_num + 1;

      cout << "  " << setw(2)  << solution_num
           << "  " << setw(10) << i << "   ";
      for ( j = 0; j < n; j++ )
      {
        cout << " " << setw(1) << bvec[j];
      }
      cout << "\n";
    }
  }
//
//  Report.
//
  cout << "\n";
  cout << "  Number of solutions found was " << solution_num << "\n";
//
//  Shut down.
//
  cout << "\n";
  cout << "SATISFY\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
# undef N
}
//****************************************************************************80

int circuit_value ( int n, int bvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCUIT_VALUE returns the value of a circuit for a given input set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Michael Quinn,
//    Parallel Programming in C with MPI and OpenMP,
//    McGraw-Hill, 2004,
//    ISBN13: 978-0071232654,
//    LC: QA76.73.C15.Q55.
//
//  Parameters:
//
//    Input, int N, the length of the input vector.
//
//    Input, int BVEC[N], the binary inputs.
//
//    Output, int CIRCUIT_VALUE, the output of the circuit.
//
{
  int value;

  value = 
       (  bvec[0]  ||  bvec[1]  )
    && ( !bvec[1]  || !bvec[3]  )
    && (  bvec[2]  ||  bvec[3]  )
    && ( !bvec[3]  || !bvec[4]  )
    && (  bvec[4]  || !bvec[5]  )
    && (  bvec[5]  || !bvec[6]  )
    && (  bvec[5]  ||  bvec[6]  )
    && (  bvec[6]  || !bvec[15] )
    && (  bvec[7]  || !bvec[8]  )
    && ( !bvec[7]  || !bvec[13] )
    && (  bvec[8]  ||  bvec[9]  )
    && (  bvec[8]  || !bvec[9]  )
    && ( !bvec[9]  || !bvec[10] )
    && (  bvec[9]  ||  bvec[11] )
    && (  bvec[10] ||  bvec[11] )
    && (  bvec[12] ||  bvec[13] )
    && (  bvec[13] || !bvec[14] )
    && (  bvec[14] ||  bvec[15] )
    && (  bvec[14] ||  bvec[16] )
    && (  bvec[17] ||  bvec[1]  )
    && (  bvec[18] || !bvec[0]  )
    && (  bvec[19] ||  bvec[1]  )
    && (  bvec[19] || !bvec[18] )
    && ( !bvec[19] || !bvec[9]  )
    && (  bvec[0]  ||  bvec[17] )
    && ( !bvec[1]  ||  bvec[20] )
    && ( !bvec[21] ||  bvec[20] )
    && ( !bvec[22] ||  bvec[20] )
    && ( !bvec[21] || !bvec[20] )
    && (  bvec[22] || !bvec[20] );

  return value;
}
//****************************************************************************80

void i4_to_bvec ( int i4, int n, int bvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_BVEC converts an integer into a binary vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, the integer.
//
//    Input, int N, the dimension of the vector.
//
//    Output, int BVEC[N], the vector of binary remainders.
//
{
  int i;

  for ( i = n - 1; 0 <= i; i-- )
  {
    bvec[i] = i4 % 2;
    i4 = i4 / 2;
  }

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
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
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
