# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>

# include "asa058.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA058_PRB.
//
//  Discussion:
//
//    ASA058_PRB tests the ASA058 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA058_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA058 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA058_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tries out the ASA058 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2008
//
//  Author:
//
//    John Burkardt
//
{
# define K 5
# define M 2
# define N 100

  int b[N];
  double d[K*M];
  double dev[K];
  double dev_sum;
  int e[K];
  int e_sum;
  double f[N];
  int i;
  ifstream input;
  char *input_filename = "points_100.txt";
  int j;
  int k2;
  int nz;
  double x[N*M];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test the CLUSTR algorithm.\n";
  cout << "  Applied Statistics Algorithm 58\n";
//
//  Read the data.
//
  cout << "\n";
  cout << "  Reading the data.\n";

  input.open ( input_filename );

  if ( !input )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  Could not open the input file: \"" << input_filename << "\"\n";
    return;
  }

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= M; j++ )
    {
      input >> x[i-1+(j-1)*N];
    }
  }

  input.close ( );
//
//  Print a few data values.
//
  cout << "\n";
  cout << "  First 5 data values:\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    cout << "  " << setw(8) << i;
    for ( j = 1; j <= M; j++ )
    {
      cout << "  " << setw(14) << x[i-1+(j-1)*N];
    }
    cout << "\n";
  }
//
//  Initialize the cluster centers arbitrarily.
//
  for ( i = 1; i <= K; i++ )
  {
    for ( j = 1; j <= M; j++ )
    {
      d[i-1+(j-1)*K] = x[i-1+(j-1)*N];
    }
  }
//
//  Compute the clusters.
//
  nz = 1;
  k2 = K;

  clustr ( x, d, dev, b, f, e, N, M, K, nz, k2 );

  cout << "\n";
  cout << "  Cluster  Population  Energy\n";
  cout << "\n";

  for ( i = 1; i <= K; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << e[i-1]
         << "  " << setw(14) << dev[i-1] << "\n";
  }

  e_sum = 0;
  dev_sum = 0.0;

  for ( i = 1; i <= K; i++ )
  {
    e_sum = e_sum + e[i-1];
    dev_sum = dev_sum + dev[i-1];
  }

  cout << "\n";
  cout << "  " << "   Total"
       << "  " << setw(8) << e_sum
       << "  " << setw(14) << dev_sum << "\n";

  return;
# undef K
# undef M
# undef N
}
