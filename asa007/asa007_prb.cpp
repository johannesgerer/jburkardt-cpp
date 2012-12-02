# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa007.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA007_PRB.
//
//  Discussion:
//
//    ASA007_PRB calls the ASA007 routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA007_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA007 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA007_PRB:\n";
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
//    TEST01 demonstrates the use of SYMINV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 15

  double a[(N_MAX*(N_MAX+1))/2];
  double afull[N_MAX*N_MAX];
  double c[(N_MAX*(N_MAX+1))/2];
  double cfull[N_MAX*N_MAX];
  double cta;
  double diff;
  int i;
  int ifault;
  int j;
  int k;
  int l;
  int n;
  int nullty;
  double w[N_MAX];

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  SYMINV computes the inverse of a positive\n";
  cout << "  definite symmetric matrix.\n";
  cout << "  A compressed storage format is used\n";
  cout << "\n";
  cout << "  Here we look at the matrix A which is\n";
  cout << "  N+1 on the diagonal and\n";
  cout << "  N   on the off diagonals.\n";

  for ( n = 1; n <= N_MAX; n++ )
  {
//
//  Set A to the lower triangle of the matrix which is N+1 on the diagonal
//  and N on the off diagonals.
//
    k = 0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j < i; j++ )
      {
        a[k] = ( double ) ( n );
        k = k + 1;
      }
      a[k] = ( double ) ( n + 1 );
      k = k + 1;
    }

    syminv ( a, n, c, w, &nullty, &ifault );

    cout << "\n";
    cout << "  Matrix order N = " << n << "\n";
    cout << "  Maxtrix nullity NULLTY = " << nullty << "\n";

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i < j; i++ )
      {
        cfull[i-1+(j-1)*n] = c[k];
        cfull[j-1+(i-1)*n] = c[k];
        k = k + 1;
      }
      cfull[j-1+(j-1)*n] = c[k];
      k = k + 1;
    }

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i < j; i++ )
      {
        afull[i-1+(j-1)*n] = a[k];
        afull[j-1+(i-1)*n] = a[k];
        k = k + 1;
      }
      afull[j-1+(j-1)*n] = a[k];
      k = k + 1;
    }
//
//  Compute C * A - I.
//
    diff = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        cta = 0.0;
        for ( k = 1; k <= n; k++ )
        {
          cta = cta + cfull[i-1+(k-1)*n] * afull[k-1+(j-1)*n];
        }
        if ( i == j )
        {
          diff = diff + pow ( 1.0 - cta, 2 );
        }
        else
        {
          diff = diff + cta * cta;
        }
      }
    }

    diff = sqrt ( diff );

    cout << "  RMS ( C * A - I ) = " << diff << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 demonstrates the use of SYMINV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 15

  double a[(N_MAX*(N_MAX+1))/2];
  double afull[N_MAX*N_MAX];
  double c[(N_MAX*(N_MAX+1))/2];
  double cfull[N_MAX*N_MAX];
  double cta;
  double diff;
  int i;
  int ifault;
  int j;
  int k;
  int l;
  int n;
  int nullty;
  double w[N_MAX];

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  SYMINV computes the inverse of a positive\n";
  cout << "  definite symmetric matrix.\n";
  cout << "  A compressed storage format is used\n";
  cout << "\n";
  cout << "  Here we look at the Hilbert matrix\n";
  cout << "  A(I,J) = 1/(I+J-1)\n";
  cout << "\n";
  cout << "  For this matrix, we expect errors to grow quickly.\n";

  for ( n = 1; n <= N_MAX; n++ )
  {
//
//  Set A to the Hilbert matrix.
//
    k = 0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        a[k] = 1.0 / ( double ) ( i + j - 1 );
        k = k + 1;
      }
    }

    syminv ( a, n, c, w, &nullty, &ifault );

    cout << "\n";
    cout << "  Matrix order N = " << n << "\n";
    cout << "  Maxtrix nullity NULLTY = " << nullty << "\n";

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i < j; i++ )
      {
        cfull[i-1+(j-1)*n] = c[k];
        cfull[j-1+(i-1)*n] = c[k];
        k = k + 1;
      }
      cfull[j-1+(j-1)*n] = c[k];
      k = k + 1;
    }

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i < j; i++ )
      {
        afull[i-1+(j-1)*n] = a[k];
        afull[j-1+(i-1)*n] = a[k];
        k = k + 1;
      }
      afull[j-1+(j-1)*n] = a[k];
      k = k + 1;
    }
//
//  Compute C * A - I.
//
    diff = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        cta = 0.0;
        for ( k = 1; k <= n; k++ )
        {
          cta = cta + cfull[i-1+(k-1)*n] * afull[k-1+(j-1)*n];
        }
        if ( i == j )
        {
          diff = diff + pow ( 1.0 - cta, 2 );
        }
        else
        {
          diff = diff + cta * cta;
        }
      }
    }

    diff = sqrt ( diff );

    cout << "  RMS ( C * A - I ) = " << diff << "\n";
  }

  return;
# undef N_MAX
}
