# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <fstream>

# include "asa136.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA136_PRB.
//
//  Discussion:
//
//    ASA136_PRB calls the ASA136 routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA136_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA136 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA136_PRB:\n";
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
//    TEST01 tries out the ASA136 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *c;
  int i;
  int *ic1;
  int ifault;
  ifstream input;
  int iter;
  int j;
  int k = 5;
  int m = 100;
  int n = 2;
  int *nc;
  int nc_sum;
  double *wss;
  double wss_sum;

  a = new double[m*n];
  c = new double[k*n];
  ic1 = new int[m];
  nc = new int[k];
  wss = new double[k];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test the KMNS algorithm,\n";
  cout << "  Applied Statistics Algorithm #136.\n";
//
//  Read the data.
//
  input.open ( "points_100.txt" );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      input >> a[i-1+(j-1)*m];
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
    for ( j = 1; j <= n; j++ )
    {
      cout << "  " << setw(14) << a[i-1+(j-1)*m];
    }
    cout << "\n";
  }
//
//  Initialize the cluster centers.
//  Here, we arbitrarily make the first K data points cluster centers.
//
  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c[i-1+(j-1)*k] = a[i-1+(j-1)*m];
    }
  }

  iter = 50;
//
//  Compute the clusters.
//
  kmns ( a, m, n, c, k, ic1, nc, iter, wss, &ifault );

  if ( ifault != 0 )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  KMNS returned IFAULT = " << ifault << "\n";
    return;
  }

  cout << "\n";
  cout << "  Cluster  Population  Energy\n";
  cout << "\n";

  nc_sum = 0;
  wss_sum = 0.0;

  for ( i = 1; i <= k; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << nc[i-1]
         << "  " << setw(14) << wss[i-1] << "\n";
    nc_sum = nc_sum + nc[i-1];
    wss_sum = wss_sum + wss[i-1];
  }

  cout << "\n";
  cout << "     Total"
       << "  " << setw(8) << nc_sum
       << "  " << setw(14) << wss_sum << "\n";

  delete [] a;
  delete [] c;
  delete [] ic1;
  delete [] nc;
  delete [] wss;

  return;
}
