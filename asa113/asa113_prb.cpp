# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <fstream>

# include "asa113.hpp"

using namespace std;

int main ( );
void test01 ( );
double crswap ( double varval[], int klass[], int clsize[], int in, int ik,
  int iv, double *critvl, int i, int j, int l, int m, int iswitch );
double crtran ( double varval[], int klass[], int clsize[], int in, int ik,
  int iv, double *critvl, int i, int m, int l, int iswitch );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA113_PRB.
//
//  Discussion:
//
//    ASA136_PRB tests the ASA113 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA113_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA113 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA113_PRB:\n";
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
//    TEST01 tries out the ASA113 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int *c;
  double *c_center;
  int *c_size;
  int ci;
  double critvl;
  int i;
  int ifault;
  ifstream input;
  int j;
  int k = 5;
  int m = 100;
  int n = 2;
  int ntrans1;
  int ntrans2;
  double *wss;

  a = new double[m*n];
  c = new int[m];
  c_center = new double[k*n];
  c_size = new int[k];
  wss = new double[k];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test the ASA113 classification algorithm.\n";
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
//  Assign points randomly to classes.
//
  for ( i = 1; i <= m; i++ )
  {
    c[i-1] = ( i % k ) + 1;
  }
//
//  Define the critical value as the sum of the squares of the distances
//  of the points to their cluster center.
//
  for ( i = 1; i <= k; i++ )
  {
    c_size[i-1] = 0;
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = 0.0;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    c_size[ci-1] = c_size[ci-1] + 1;
    for ( j = 1; j <= n; j++ )
    {
      c_center[ci-1+(j-1)*k] = c_center[ci-1+(j-1)*k] + a[i-1+(j-1)*m];
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = c_center[i-1+(j-1)*k] / ( double ) ( c_size[i-1] );
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    wss[i-1] = 0.0;
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    for ( j = 1; j <= n; j++ )
    {
      wss[ci-1] = wss[ci-1] + pow ( a[i-1+(j-1)*m] - c_center[ci-1+(j-1)*k], 2 );
    }
  }

  critvl = 0.0;
  for ( i = 1; i <= k; i++ )
  {
    critvl = critvl + wss[i-1];
  }

  cout << "\n";
  cout << "        Initial CRITVL = " << critvl << "\n";
//
//  Compute the clusters.
//
  ntrans1 = -1;
  ntrans2 = -1;

  for ( ; ; )
  {
    trnsfr ( a, c, c_size, m, k, n, &critvl, &ntrans1, &ifault );

    if ( ntrans1 == 0 && ntrans2 == 0 )
    {
      break;
    }

    cout << "  After TRNSFR, CRITVL = " << critvl << "\n";

    swap ( a, c, c_size, m, k, n, &critvl, &ntrans2, &ifault );

    if ( ntrans1 == 0 && ntrans2 == 0 )
    {
      break;
    }

    cout << "    After SWAP, CRITVL = " << critvl << "\n";
  }
//
//  Define the critical value as the sum of the squares of the distances
//  of the points to their cluster center.
//
  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = 0.0;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    for ( j = 1; j <= n; j++ )
    {
      c_center[ci-1+(j-1)*k] = c_center[ci-1+(j-1)*k] + a[i-1+(j-1)*m];
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = c_center[i-1+(j-1)*k] / ( double ) ( c_size[i-1] );
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    wss[i-1] = 0.0;
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    for ( j = 1; j <= n; j++ )
    {
      wss[ci-1] = wss[ci-1]
                + pow ( a[i-1+(j-1)*m] - c_center[ci-1+(j-1)*k], 2 );
    }
  }
  cout << "\n";
  cout << "  Cluster  Population  Energy\n";
  cout << "\n";

  for ( i = 1; i <= k; i++ )
  {
  cout << "  " << setw(8)  << i
       << "  " << setw(8)  << c_size[i-1]
       << "  " << setw(14) << wss[i-1] << "\n";
  }
  cout << "\n";
  cout << "     Total"
       << "  " << setw(8)  << m
       << "  " << setw(14) << critvl << "\n";

  delete [] a;
  delete [] c;
  delete [] c_center;
  delete [] c_size;
  delete [] wss;

  return;
}
//****************************************************************************80

double crswap ( double a[], int c[], int c_size[], int m, int k,
  int n, double *critvl, int i1, int i2, int c1, int c2, int iswitch )

//****************************************************************************80
//
//  Purpose:
//
//    CRSWAP determines the effect of swapping two objects.
//
//  Discussion:
//
//    This computation is very inefficient.  It is only set up so that we
//    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
//    ASA 136.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Colin Banfield, LC Bassill,
//    Algorithm AS 113:
//    A transfer for non-hierarchichal classification,
//    Applied Statistics,
//    Volume 26, Number 2, 1977, pages 206-210.
//
//  Parameters:
//
//    Input, double A(M,N), the data values.  There are M objects,
//    each having spatial dimension N.
//
//    Input, int C(M), the classification of each object.
//
//    Input, int C_SIZE(K), the number of objects in each class.
//
//    Input, int M, the number of objects.
//
//    Input, int K, the number of classes.
//
//    Input, int N, the number of spatial dimensions, or variates,
//    of the objects.
//
//    Input, double *CRITVL, the current value of the criterion.
//
//    Input, int I1, I2, the objects to be swapped.
//
//    Input, int C1, C2, the current classes of objects I1 and I2.
//
//    Input, int ISWITCH:
//    1, indicates that I1 and I2 should be temporarily swapped, the
//       change in CRITVL should be computed, and then I1 and I2 restored.
//    2, indicates that I1 and I2 will be swapped.
//
//    Output, double CRSWAP, the change to CRITVL that would occur if I1 and
//    I2 were swapped.  This is only computed for ISWITCH = 1.
//
{
  double *c_center;
  int ci;
  double critvl_new;
  int i;
  double inc;
  int j;

  if ( iswitch == 2 )
  {
    inc = 0.0;
    return inc;
  }

  c_center = new double[k*n];
//
//  Move object I1 from class C1 to class C2.
//  Move object I2 from class C2 to class C1.
//
  c[i1-1] = c2;
  c[i2-1] = c1;
//
//  Define the critical value as the sum of the squares of the distances
//  of the points to their cluster center.
//
  for ( i = 1; i <= k; i++ )
  {
    c_size[i-1] = 0;
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = 0.0;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    c_size[ci-1] = c_size[ci-1] + 1;
    for ( j = 1; j <= n; j++ )
    {
      c_center[ci-1+(j-1)*k] = c_center[ci-1+(j-1)*k] + a[i-1+(j-1)*m];
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = c_center[i-1+(j-1)*k] / ( double ) ( c_size[i-1] );
    }
  }

  critvl_new = 0.0;

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    for ( j = 1; j <= n; j++ )
    {
      critvl_new = critvl_new
                 + pow ( a[i-1+(j-1)*m] - c_center[ci-1+(j-1)*k], 2 );
    }
  }

  inc = critvl_new - *critvl;
//
//  Move object I1 from class C2 to class C1.
//  Move object I2 from class C1 to class C2.
//
  c[i1-1] = c1;
  c[i2-1] = c2;

  delete [] c_center;

  return inc;
}
//****************************************************************************80

double crtran ( double a[], int c[], int c_size[], int m, int k, int n,
  double *critvl, int i1, int c1, int c2, int iswitch )

//****************************************************************************80
//
//  Purpose:
//
//    CRTRAN determines the effect of moving an object to another class.
//
//  Discussion:
//
//    This computation is very inefficient.  It is only set up so that we
//    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
//    ASA 136.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Colin Banfield, LC Bassill,
//    Algorithm AS 113:
//    A transfer for non-hierarchichal classification,
//    Applied Statistics,
//    Volume 26, Number 2, 1977, pages 206-210.
//
//  Parameters:
//
//    Input, double AL(M,N), the data values.  There are M objects,
//    each having spatial dimension N.
//
//    Input, int C(M), the classification of each object.
//
//    Input, int C_SIZE(K), the number of objects in each class.
//
//    Input, int M, the number of objects.
//
//    Input, int K, the number of classes.
//
//    Input, int N, the number of spatial dimensions, or variates,
//    of the objects.
//
//    Input, double *CRITVL, the current value of the criterion.
//
//    Input, int I1, the object to be transferred.
//
//    Input, int C1, C2, the current class of object I1, and the
//    class to which it may be transferred.
//
//    Input, int ISWITCH:
//    1, indicates that I1 should be temporarily transferred, the
//       change in CRITVL should be computed, and then I1 restored.
//    2, indicates that I1 will be permanently transferred.
//
//    Output, double CRTRAN, the change to CRITVL that would occur if I1 were
//    transferred from class C1 to C2.  This is only computed for ISWITCH = 1.
//
{
  double *c_center;
  int ci;
  double critvl_new;
  int i;
  double inc;
  int j;

  if ( iswitch == 2 )
  {
    inc = 0.0;
    return inc;
  }

  c_center = new double[k*n];
//
//  Move object I from class C1 to class C2.
//
  c[i1-1] = c2;
  c_size[c1-1] = c_size[c1-1] - 1;
  c_size[c2-1] = c_size[c2-1] + 1;
//
//  Define the critical value as the sum of the squares of the distances
//  of the points to their cluster center.
//
  for ( i = 1; i <= k; i++ )
  {
    c_size[i-1] = 0;
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = 0.0;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    c_size[ci-1] = c_size[ci-1] + 1;
    for ( j = 1; j <= n; j++ )
    {
      c_center[ci-1+(j-1)*k] = c_center[ci-1+(j-1)*k] + a[i-1+(j-1)*m];
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = c_center[i-1+(j-1)*k] / ( double ) ( c_size[i-1] );
    }
  }

  critvl_new = 0.0;

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    for ( j = 1; j <= n; j++ )
    {
      critvl_new = critvl_new
                 + pow ( a[i-1+(j-1)*m] - c_center[ci-1+(j-1)*k], 2 );
    }
  }

  inc = critvl_new - *critvl;
//
//  Move object I1 from class C2 to class C1.
//
  c[i1-1] = c1;
  c_size[c1-1] = c_size[c1-1] + 1;
  c_size[c2-1] = c_size[c2-1] - 1;

  delete [] c_center;

  return inc;
}
