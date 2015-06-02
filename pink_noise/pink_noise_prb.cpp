# include <cstdlib>
# include <iostream>
# include <iomanip>

# include "pink_noise.hpp"

using namespace std;

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
//    MAIN is the main program for PINK_NOISE_PRB.
//
//  Discussion:
//
//    PINK_NOISE_PRB tests the PINK_NOISE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "PINK_NOISE_PRB:\n";
  cout << "  Test the PINK_NOISE library.\n";
  cout << "  C++ version\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PINK_NOISE_PRB:\n";
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
//    TEST01 tests WRAP2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int m;
  int q;
  int q_in;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  WRAP2 performs a circular wrap.\n";
  cout << "  Q is expected to range between 0 and M.\n";
  cout << "  WRAP2 takes an input value of Q, and either\n";
  cout << "  increments it by M+1 until in the range, or\n";
  cout << "  decrements it by M+1 until in the range,\n";
  cout << "  and returns the result as the function value.\n";

  for ( m = 2; m <= 4; m++ )
  {
    cout << "\n";
    cout << "   M  Qin  Qout\n";
    cout << "\n";
    for ( i = -5; i < 3 * m; i++ )
    {
      q = i;
      q_in = q;
      wrap2 ( m, &q );
      cout << "  " << setw(2) << m
           << "  " << setw(2) << q_in
           << "  " << setw(2) << q << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CDELAY2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int m;
  int q;
  int q_in;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  CDELAY2 is a circular buffer implementation\n";
  cout << "  of an M-fold delay.  Q is a counter\n";
  cout << "  which is decremented by CDELAY2, but reset to M\n";
  cout << "  after it reaches 0.\n";

  for ( m = 2; m <= 4; m++ )
  {
    cout << "\n";
    cout << "   I   M  Qin  Qout\n";
    cout << "\n";
    q = m;
    for ( i = 1; i <= 3 * ( m + 1 ); i++ )
    {
      q_in = q;
      cdelay2 ( m, &q );
      cout << "  " << setw(2) << i
           << "  " << setw(2) << m
           << "  " << setw(2) << q_in
           << "  " << setw(2) << q << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests RANH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i;
  int q;
  double u;
  double y;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  RANH is a random hold function.\n";
  cout << "  Given a value U and a delay D, it returns the value\n";
  cout << "  U for D calls, then resets U.\n";

  for ( d = 5; 1 <= d; d-- )
  {
    cout << "\n";
    cout << "   I   D   Q      U           Y\n";
    cout << "\n";
    u = 0.5;
    q = 3;
    for ( i = 1; i <= 20; i++ )
    {
      y = ranh ( d, &u, &q );
      cout << "  " << setw(2) << i
           << "  " << setw(2) << d
           << "  " << setw(2) << q 
           << "  " << setw(10) << u
           << "  " << setw(10) << y << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests RAN1F.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  int b;
  int i;
  int *q;
  int rep;
  double *u;
  double y;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  RAN1F generates random values with an approximate\n";
  cout << "  1/F distribution.\n";

  for ( b = 1; b < 32; b = b * 2 )
  {
    u = new double[b];
    q = new int[b];
    for ( rep = 1; rep <= 4; rep++ )
    {
      for ( i = 0; i < b; i++ )
      {
        u[i] = ( double ) rand ( ) / ( double ) ( RAND_MAX ) - 0.5;
      }
      for ( i = 0; i < b; i++ )
      {
        q[i] = 0;
      }
      cout << "\n";
      cout << "   B   I      Y\n";
      cout << "\n";

      for ( i = 1; i <= 20; i++ )
      {
        y = ran1f ( b, u, q );
        cout << "  " << setw(2) << b
             << "  " << setw(2) << i
             << "  " << setw(10) << y << "\n";
      }
    }
    delete [] q;
    delete [] u;
  }
  return;
}
