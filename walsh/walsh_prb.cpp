# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "walsh.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    WALSH_PRB calls the WALSH test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "WALSH_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the WALSH library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "WALSH_PRB\n";
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
//    TEST01 tests FWT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n = 16;
  int seed;
  double *w;
  double *x;
  double *y;
  double *z;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  FWT computes a fast Walsh transform.\n";

  for ( j = 1; j <= 2; j++ )
  {
    if ( j == 1 )
    {
      seed = 123456789;
      w = r8vec_uniform_01_new ( n, &seed );
    }
    else
    {
      w = new double[n];
      for ( i = 0; i < n; i++ )
      {
        w[i] = ( double ) ( i + 1 );
      }
    }

    x = r8vec_copy_new ( n, w );
    fwt ( n, w );
    y = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      y[i] = y[i] / ( double ) ( n );
    }
    fwt ( n, w );
    z = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      z[i] = z[i] / ( double ) ( n );
    }

    cout << "\n";
    cout << "     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N\n";
    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << i
           << "  " << setw(10) << x[i]
           << "  " << setw(10) << y[i]
           << "  " << setw(10) << z[i] << "\n";
    }
    delete [] w;
    delete [] x;
    delete [] y;
    delete [] z;
  }
    
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests WALSH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n = 16;
  int seed;
  double *w;
  double *x;
  double *y;
  double *z;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  WALSH computes a fast Walsh transform.\n";

  for ( j = 1; j <= 2; j++ )
  {
    if ( j == 1 )
    {
       seed = 123456789;
       w = r8vec_uniform_01_new ( n, &seed );
    }
    else
    {
      w = new double[n];
      for ( i = 0; i < n; i++ )
      {
        w[i] = ( double ) ( i + 1 );
      }
    }

    x = r8vec_copy_new ( n, w );
    walsh ( n, w );
    y = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      y[i] = y[i] / ( double ) ( n );
    }
    walsh ( n, w );
    z = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      z[i] = z[i] / ( double ) ( n );
    }

    cout << "\n";
    cout << "     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N\n";
    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << i
           << "  " << setw(10) << x[i]
           << "  " << setw(10) << y[i]
           << "  " << setw(10) << z[i] << "\n";
    }

    delete [] w;
    delete [] x;
    delete [] y;
    delete [] z;
  }
    
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests HAAR, HAARIN and HNORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n = 16;
  int seed;
  double *w;
  double *x;
  double *y;
  double *z;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  HAAR computes a Haar transform.\n";
  cout << "  HNORM normalizes the transformed data.\n";
  cout << "  HAARIN computes an inverse Haar transform.\n";

  for ( j = 1; j <= 2; j++ )
  {
    if ( j == 1 )
    {
      seed = 123456789;
      w = r8vec_uniform_01_new ( n, &seed );
    }
    else
    {
      w = new double[n];
      for ( i = 0; i < n; i++ )
      {
        w[i] = ( double ) ( i + 1 );
      }
    }

    x = r8vec_copy_new ( n, w );

    haar ( n, w );

    y = r8vec_copy_new ( n, w );

    hnorm ( n, w );

    z = r8vec_copy_new ( n, w );

    haarin ( n, w );

    cout << "\n";
    cout << "     I        X(I)    Y=HAAR(X)  Z=HNORM(Y)  W=HAARIN(Z)\n";
    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << i
           << "  " << setw(10) << x[i]
           << "  " << setw(10) << y[i]
           << "  " << setw(10) << z[i]
           << "  " << setw(10) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
    delete [] y;
    delete [] z;
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests FFWT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n = 16;
  int seed;
  double *w;
  double *x;
  double *y;
  double *z;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  FFWT computes a fast Walsh transform.\n";

  for ( j = 1; j <= 2; j++ )
  {
    if ( j == 1 )
    {
      seed = 123456789;
      w = r8vec_uniform_01_new ( n, &seed );
    }
    else
    {
      w = new double[n];
      for ( i = 0; i < n; i++ )
      {
        w[i] = ( double ) ( i + 1 );
      }
    }
    x = r8vec_copy_new ( n, w );
    ffwt ( n, w );
    y = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      y[i] = y[i] / ( double ) ( n );
    }
    ffwt ( n, w );
    z = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      z[i] = z[i] / ( double ) ( n );
    }

    cout << "\n";
    cout << "     I        X(I)   Y=FFWT(X)/N  Z=FFWT(Y)/N\n";
    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(4) << i
           << "  " << setw(10) << x[i]
           << "  " << setw(10) << y[i]
           << "  " << setw(10) << z[i] << "\n";
    }
    delete [] w;
    delete [] x;
    delete [] y;
    delete [] z;
  } 
  return;
}
