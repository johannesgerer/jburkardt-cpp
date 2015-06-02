# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa006.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA006_PRB.
//
//  Discussion:
//
//    ASA006_PRB tests the ASA006 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA006_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA006 library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA006_PRB:\n";
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
//    TEST01 demonstrates the use of CHOLESKY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 15

  double a[(N_MAX*(N_MAX+1))/2];
  double diff;
  int i;
  int ifault;
  int j;
  int k;
  int l;
  int n;
  int nn;
  int nullty;
  double r[N_MAX];
  double rmax;
  double u[(N_MAX*(N_MAX+1))/2];
  double ufull[N_MAX*N_MAX];
  double utu;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  CHOLESKY computes the Cholesky factorization\n";
  cout << "  of a positive definite symmetric matrix.\n";
  cout << "  A compressed storage format is used\n";
  cout << "\n";
  cout << "  Here we look at the matrix A which is\n";
  cout << "  N+1 on the diagonal and\n";
  cout << "  N   on the off diagonals.\n";

  for ( n = 1; n <= N_MAX; n++ )
  {
    nn = ( n * ( n + 1 ) ) / 2;
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

    cholesky ( a, n, nn, u, &nullty, &ifault );

    cout << "\n";
    cout << "  Matrix order N = " << n << "\n";
    cout << "  Maxtrix nullity NULLTY = " << nullty << "\n";

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= j; i++ )
      {
        ufull[i-1+(j-1)*n] = u[k];
        k = k + 1;
      }
      for ( i = j + 1; i <= n; i++ )
      {
        ufull[i-1+(j-1)*n] = 0.0;
      }
    }
//
//  Compute U' * U and compare to A.
//
    k = 0;
    diff = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        utu = 0.0;
        for ( l = 1; l <= n; l++ )
        {
          utu = utu + ufull[l-1+(i-1)*n] * ufull[l-1+(j-1)*n];
        }
        diff = diff + ( a[k] - utu ) * ( a[k] - utu );
        k = k + 1;
      }
    }

    diff = sqrt ( diff );

    cout << "  RMS ( A - U'*U ) = " << diff << "\n";
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
//    TEST02 demonstrates the use of CHOLESKY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 15

  double a[(N_MAX*(N_MAX+1))/2];
  double diff;
  int i;
  int ifault;
  int j;
  int k;
  int l;
  int n;
  int nn;
  int nullty;
  double u[(N_MAX*(N_MAX+1))/2];
  double ufull[N_MAX*N_MAX];
  double utu;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  CHOLESKY computes the Cholesky factorization\n";
  cout << "  of a positive definite symmetric matrix.\n";
  cout << "  A compressed storage format is used\n";
  cout << "\n";
  cout << "  Here we look at the Hilbert matrix\n";
  cout << "  A(I,J) = 1/(I+J-1)\n";
  cout << "\n";
  cout << "  For this matrix, we expect errors to grow quickly.\n";

  for ( n = 1; n <= N_MAX; n++ )
  {
    nn = ( n * ( n + 1 ) ) / 2;
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

    cholesky ( a, n, nn, u, &nullty, &ifault );

    cout << "\n";
    cout << "  Matrix order N = " << n << "\n";
    cout << "  Maxtrix nullity NULLTY = " << nullty << "\n";

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= j; i++ )
      {
        ufull[i-1+(j-1)*n] = u[k];
        k = k + 1;
      }
      for ( i = j + 1; i <= n; i++ )
      {
        ufull[i-1+(j-1)*n] = 0.0;
      }
    }
//
//  Compute U' * U and compare to A.
//
    k = 0;
    diff = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        utu = 0.0;
        for ( l = 1; l <= n; l++ )
        {
          utu = utu + ufull[l-1+(i-1)*n] * ufull[l-1+(j-1)*n];
        }
        diff = diff + ( a[k] - utu ) * ( a[k] - utu );
        k = k + 1;
      }
    }

    diff = sqrt ( diff );

    cout << "  RMS ( A - U'*U ) = " << diff << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 demonstrates the use of SUBCHL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 15
# define NN_MAX ((N_MAX*(N_MAX+1))/2)

  double a[NN_MAX];
  int b[N_MAX];
  double det;
  double diff;
  int i;
  int ifault;
  int j;
  int k;
  int l;
  int n;
  int nullty;
  double u[NN_MAX];
  double ufull[N_MAX*N_MAX];
  double utu;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  SUBCHL computes the Cholesky factor\n";
  cout << "  of a submatrix\n";
  cout << "  of a positive definite symmetric matrix.\n";
  cout << "  A compressed storage format is used.\n";
  cout << "\n";
  cout << "  Here we look at the Hilbert matrix\n";
  cout << "  A(I,J) = 1/(I+J-1).\n";
  cout << "\n";
  cout << "  For this particular matrix, we expect the\n";
  cout << "  errors to grow rapidly.\n";
//
//  Set A to the N_MAX order Hilbert matrix.
//
  k = 0;
  for ( i = 1; i <= N_MAX; i++ )
  {
    for ( j = 1; j <= i; j++ )
    {

      a[k] = 1.0 / ( double ) ( i + j - 1 );
      k = k + 1;
    }
  }

  for ( n = 1; n <= N_MAX; n++ )
  {
    for ( i = 1; i <= n; i++ )
    {
      b[i-1] = i;
    }

    subchl ( a, b, n, u, &nullty, &ifault, NN_MAX, &det );

    cout << "\n";
    cout << "  Matrix order N = " << n << "\n";
    cout << "  Maxtrix nullity NULLTY = " << nullty << "\n";
    cout << "  Matrix determinant DET = " << det << "\n";

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= j; i++ )
      {
        k = k + 1;
        ufull[i-1+(j-1)*n] = u[k-1];
      }
      for ( i = j + 1; i <= n; i++ )
      {
        ufull[i-1+(j-1)*n] = 0.0;
      }
    }
//
//  Compute U' * U and compare to A.
//
    k = 0;
    diff = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        k = k + 1;
        utu = 0.0;
        for ( l = 1; l <= n; l++ )
        {
          utu = utu + ufull[l-1+(i-1)*n] * ufull[l-1+(j-1)*n];
        }
        diff = diff + pow ( a[k-1] - utu, 2 );
      }
    }
    diff = sqrt ( diff );
    cout << "  RMS ( A - U'*U ) = " << diff << "\n";
  }

  return;
# undef N_MAX
# undef NN_MAX
}
