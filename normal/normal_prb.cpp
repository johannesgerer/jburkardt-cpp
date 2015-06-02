# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <complex>

using namespace std;

# include "normal.hpp"

int main ( );
void c4_normal_01_test ( );
void c8_normal_01_test ( );
void i4_normal_ab_test ( );
void i8_normal_ab_test ( );
void r4_normal_01_test ( );
void r4_normal_ab_test ( );
void r8_normal_01_test ( );
void r8_normal_ab_test ( );
void r8mat_normal_01_new_test ( );
void r8vec_normal_01_new_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for NORMAL_PRB.
//
//  Discussion:
//
//    NORMAL_PRB tests the NORMAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "NORMAL_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the NORMAL library\n";

  c4_normal_01_test ( );
  c8_normal_01_test ( );
  i4_normal_ab_test ( );
  i8_normal_ab_test ( );
  r4_normal_01_test ( );
  r4_normal_ab_test ( );
  r8_normal_01_test ( );
  r8_normal_ab_test ( );
  r8mat_normal_01_new_test ( );
  r8vec_normal_01_new_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "NORMAL_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void c4_normal_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NORMAL_01_TEST tests C4_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  complex <float> r;
  int seed;

  cout << "\n";
  cout << "C4_NORMAL_01_TEST\n";
  cout << "  C4_NORMAL_01 computes pseudorandom complex values\n";
  cout << "  normally distributed in the unit circle.\n";

  seed = 123456789;

  cout << "\n";
  cout << "  SEED = " << seed << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    r = c4_normal_01 ( seed );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << real ( r )
         << "  " << setw(12) << imag ( r ) << "\n";
  }

  return;
}
//****************************************************************************80

void c8_normal_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_NORMAL_01_TEST tests C8_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  complex <double> r;
  int seed;

  cout << "\n";
  cout << "C8_NORMAL_01_TEST\n";
  cout << "  C8_NORMAL_01 computes pseudorandom complex values\n";
  cout << "  normally distributed in the unit circle.\n";

  seed = 123456789;

  cout << "\n";
  cout << "  SEED = " << seed << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    r = c8_normal_01 ( seed );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << real ( r )
         << "  " << setw(12) << imag ( r ) << "\n";
  }

  return;
}
//****************************************************************************80

void i4_normal_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_NORMAL_AB_TEST tests I4_NORMAL_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  float mu;
  int r;
  int seed;
  float sigma;

  cout << "\n";
  cout << "I4_NORMAL_AB_TEST\n";
  cout << "  I4_NORMAL_AB computes pseudonormal integer values\n";
  cout << "  with mean MU and standard deviation SIGMA.\n";

  mu = 70.0;
  sigma = 10.0;
  seed = 123456789;

  cout << "\n";
  cout << "  The mean = " << mu << "\n";
  cout << "  The standard deviation = " << sigma << "\n";
  cout << "  SEED = " << seed << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    r = i4_normal_ab ( mu, sigma, seed );
    cout << "  " << setw(8) << i
         << "  " << setw(8) << r << "\n";
  }

  return;
}
//****************************************************************************80

void i8_normal_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I8_NORMAL_AB_TEST tests I8_NORMAL_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double mu;
  long long int r;
  long long int seed;
  double sigma;

  cout << "\n";
  cout << "I8_NORMAL_AB_TEST\n";
  cout << "  I8_NORMAL_AB computes pseudonormal integer values\n";
  cout << "  with mean MU and standard deviation SIGMA.\n";

  mu = 70.0;
  sigma = 10.0;
  seed = 123456789;

  cout << "\n";
  cout << "  The mean = " << mu << "\n";
  cout << "  The standard deviation = " << sigma << "\n";
  cout << "  SEED = " << seed << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    r = i8_normal_ab ( mu, sigma, seed );
    cout << "  " << setw(8) << i
         << "  " << setw(8) << r << "\n";
  }

  return;
}
//****************************************************************************80

void r4_normal_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NORMAL_01_TEST tests R4_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  float r;
  int seed;
 
  cout << "\n";
  cout << "R4_NORMAL_01_TEST\n";
  cout << "  R4_NORMAL_01 computes pseudonormal values\n";
  cout << "  with mean 0.0 and standard deviation 1.0.\n";

  seed = 123456789;

  cout << "\n";
  cout << "  SEED = " << seed << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    r = r4_normal_01 ( seed );
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r << "\n";
  }

  return;
}
//****************************************************************************80

void r4_normal_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NORMAL_AB_TEST tests R4_NORMAL_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  float mu;
  float r;
  int seed;
  float sigma;

  cout << "\n";
  cout << "R4_NORMAL_AB_TEST\n";
  cout << "  R4_NORMAL_AB computes pseudonormal values\n";
  cout << "  with mean MU and standard deviation SIGMA.\n";

  mu = 10.0;
  sigma = 2.0;
  seed = 123456789;

  cout << "\n";
  cout << "  The mean = " << mu << "\n";
  cout << "  The standard deviation = " << sigma << "\n";
  cout << "  SEED = " << seed << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    r = r4_normal_ab ( mu, sigma, seed );
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r << "\n";
  }

  return;
}
//****************************************************************************80

void r8_normal_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01_TEST tests R8_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double r;
  int seed;

  cout << "\n";
  cout << "R8_NORMAL_01_TEST\n";
  cout << "  R8_NORMAL_01 computes pseudonormal values\n";
  cout << "  with mean 0.0 and standard deviation 1.0.\n";

  seed = 123456789;

  cout << "\n";
  cout << "  SEED = " << seed << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    r = r8_normal_01 ( seed );
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r << "\n";
  }

  return;
}
//****************************************************************************80

void r8_normal_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_AB_TEST tests R8_NORMAL_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double mu;
  double r;
  int seed;
  double sigma;

  cout << "\n";
  cout << "R8_NORMAL_AB_TEST\n";
  cout << "  R8_NORMAL_AB computes pseudonormal values\n";
  cout << "  with mean MU and standard deviation SIGMA.\n";

  mu = 10.0;
  sigma = 2.0;
  seed = 123456789;

  cout << "\n";
  cout << "  The mean = " << mu << "\n";
  cout << "  The standard deviation = " << sigma << "\n";
  cout << "  SEED = " << seed << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    r = r8_normal_ab ( mu, sigma, seed );
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r << "\n";
  }

  return;
}
//****************************************************************************80

void r8mat_normal_01_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORMAL_01_NEW_TEST tests R8MAT_NORMAL_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  double *r;
  int seed;

  cout << "\n";
  cout << "R8MAT_NORMAL_01_NEW_TEST\n";
  cout << "  R8MAT_NORMAL_01_NEW computes a matrix of values.\n";

  m = 5;
  n = 4;
  seed = 123456789;
  cout << "\n";
  cout << "  SEED = " << seed << "\n";

  r = r8mat_normal_01_new ( m, n, seed );

  r8mat_print ( m, n, r, "  Matrix:" );
  
  delete [] r;

  return;
}
//****************************************************************************80

void r8vec_normal_01_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_NEW_TEST tests R8VEC_NORMAL_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int n = 10;
  double *r;
  int seed;

  cout << "\n";
  cout << "R8VEC_NORMAL_01_NEW_TEST\n";
  cout << "  R8VEC_NORMAL_01_NEW computes a vector of Normal 01 values.\n";

  seed = 123456789;

  cout << "\n";
  cout << "  SEED = " << seed << "\n";

  r = r8vec_normal_01_new ( n, seed );

  r8vec_print ( n, r, "  Random vector:" );
  
  delete [] r;

  return;
}


