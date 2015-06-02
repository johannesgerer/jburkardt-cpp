# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <string>

using namespace std;

# include "polynomial.hpp"

int main ( );
void polynomial_add_test ( );
void polynomial_axpy_test ( );
void polynomial_compress_test ( );
void polynomial_dif_test ( );
void polynomial_mul_test ( );
void polynomial_print_test ( );
void polynomial_scale_test ( );
void polynomial_sort_test ( );
void polynomial_value_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for POLYNOMIAL_PRB.
//
//  Discussion:
//
//    POLYNOMIAL_PRB tests the POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "POLYNOMIAL_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the POLYNOMIAL library.\n";

  polynomial_add_test ( );
  polynomial_axpy_test ( );
  polynomial_compress_test ( );
  polynomial_dif_test ( );
  polynomial_mul_test ( );
  polynomial_print_test ( );
  polynomial_scale_test ( );
  polynomial_sort_test ( );
  polynomial_value_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "POLYNOMIAL_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void polynomial_add_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_ADD_TEST tests POLYNOMIAL_ADD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  double c[11];
  double c1[6] = { 7.0, - 5.0, 9.0, 11.0, 0.0, - 13.0 };
  double c2[5] = { 2.0, 3.0, -8.0, 4.0, 9.0 };
  int m = 3;
  int e[11];
  int e1[6] = { 1, 2, 4, 5, 12, 33 };
  int e2[5] = { 1, 3, 4, 30, 33 };
  int o;
  int o1 = 6;
  int o2 = 5;
  int o_max = o1 + o2;
  string title = "  P1(X) + P2(X) =";
  string title1 = "  P1(X) =";
  string title2 = "  P2(X) =";

  cout << "\n";
  cout << "POLYNOMIAL_ADD_TEST\n";
  cout << "  POLYNOMIAL_ADD adds two polynomials.\n";

  cout << "\n";
  polynomial_print ( m, o1, c1, e1, title1 );

  cout << "\n";
  polynomial_print ( m, o2, c2, e2, title2 );

  polynomial_add ( o1, c1, e1, o2, c2, e2, o, c, e );
  cout << "\n";
  polynomial_print ( m, o, c, e, title );

  return;
}
//****************************************************************************80

void polynomial_axpy_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_AXPY_TEST tests POLYNOMIAL_ADD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c[11];
  double c1[6] = { 7.0, - 5.0, 9.0, 11.0, 0.0, - 13.0 };
  double c2[5] = { 2.0, 3.0, -8.0, 4.0, 9.0 };
  int m = 3;
  int e[11];
  int e1[6] = { 1, 2, 4, 5, 12, 33 };
  int e2[5] = { 1, 3, 4, 30, 33 };
  int o;
  int o1 = 6;
  int o2 = 5;
  int o_max = o1 + o2;
  double s = 10.0;
  string title = "  S * P1(X) + P2(X) =";
  string title1 = "  P1(X) =";
  string title2 = "  P2(X) =";

  cout << "\n";
  cout << "POLYNOMIAL_AXPY_TEST\n";
  cout << "  POLYNOMIAL_AXPY adds a multiple of one polynomial to another.\n";

  cout << "\n";
  polynomial_print ( m, o1, c1, e1, title1 );

  cout << "\n";
  polynomial_print ( m, o2, c2, e2, title2 );

  polynomial_axpy ( s, o1, c1, e1, o2, c2, e2, o, c, e );
  cout << "\n";
  polynomial_print ( m, o, c, e, title );

  return;
}
//****************************************************************************80

void polynomial_compress_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_COMPRESS_TEST tests POLYNOMIAL_COMPRESS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c[10] = { 7.0, - 5.0, 5.0, 9.0, 11.0, 3.0, 6.0, 0.0, - 13.0, 1.0E-20 };
  double c2[10];
  int m = 3;
  int e[10] = { 1, 2, 2, 4, 5, 5, 5, 12, 33, 35 }; 
  int e2[10];
  int j;
  int nx = 2;
  int o = 10;
  int o2;
  string title;

  cout << "\n";
  cout << "POLYNOMIAL_COMPRESS_TEST\n";
  cout << "  POLYNOMIAL_COMPRESS compresses a polynomial.\n";

  cout << "\n";
  title = "  Uncompressed P(X) = ";
  polynomial_print ( m, o, c, e, title );

  polynomial_compress ( o, c, e, o2, c2, e2 );

  cout << "\n";
  title = "  Compressed P(X) = ";
  polynomial_print ( m, o2, c2, e2, title );

  return;
}
//****************************************************************************80

void polynomial_dif_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_DIF_TEST tests POLYNOMIAL_DIF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  double c[4];
  double c1[4] = { 2.0, 3.0, 4.0, 5.0 };
  int m = 2;
  int dif[2] = { 2, 1 };
  int e[4];
  int e1[4] = { 1, 10, 12, 32 };
  int o;
  int o1 = 4;
  int o_max = 4;
  string title = "  d3 P(X) dx1 dx1 dx2 =";
  string title1 = "  P(X) =";

  cout << "\n";
  cout << "POLYNOMIAL_DIF_TEST\n";
  cout << "  POLYNOMIAL_DIF computes derivatives of a polynomial.\n";

  cout << "\n";
  polynomial_print ( m, o1, c1, e1, title1 );

  polynomial_dif ( m, o1, c1, e1, dif, o, c, e );
  cout << "\n";
  polynomial_print ( m, o, c, e, title );

  return;
}
//****************************************************************************80

void polynomial_mul_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_MUL_TEST tests POLYNOMIAL_MUL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  double c[8];
  double c1[4] = { 2.0, 3.0, 4.0, 5.0 };
  double c2[2] = { 6.0, 7.0 };
  int m = 3;
  int e[8];
  int e1[4] = { 1, 3, 4, 6 };
  int e2[2] = { 2, 5 };
  int o;
  int o1 = 4;
  int o2 = 2;
  int o_max = o1 * o2;
  string title = "  P1(X) * P2(X) =";
  string title1 = "  P1(X) =";
  string title2 = "  P2(X) =";

  cout << "\n";
  cout << "POLYNOMIAL_MUL_TEST\n";
  cout << "  POLYNOMIAL_MUL multiplies two polynomials.\n";

  cout << "\n";
  polynomial_print ( m, o1, c1, e1, title1 );

  cout << "\n";
  polynomial_print ( m, o2, c2, e2, title2 );

  polynomial_mul ( m, o1, c1, e1, o2, c2, e2, o, c, e );
  cout << "\n";
  polynomial_print ( m, o, c, e, title );

  return;
}
//****************************************************************************80

void polynomial_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_PRINT_TEST tests POLYNOMIAL_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c[6] = {7.0, - 5.0, 9.0, 11.0, 0.0, - 13.0 };
  int m = 3;
  int e[6] = { 1, 2, 4, 5, 12, 33 };
  int o = 6;
  string title = "  P1(X) =";

  cout << "\n";
  cout << "POLYNOMIAL_PRINT_TEST\n";
  cout << "  POLYNOMIAL_PRINT prints a polynomial.\n";

  cout << "\n";
  polynomial_print ( m, o, c, e, title );

  return;
}
//****************************************************************************80

void polynomial_scale_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_SCALE_TEST tests POLYNOMIAL_SCALE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c[6] = {7.0, - 5.0, 9.0, 11.0, 0.0, - 13.0 };
  int m = 3;
  int e[6] = { 1, 2, 4, 5, 12, 33 };
  int o = 6;
  double s;
  string title;

  cout << "\n";
  cout << "POLYNOMIAL_PRINT_TEST\n";
  cout << "  POLYNOMIAL_PRINT prints a polynomial.\n";

  cout << "\n";
  title = "  Original P(X):";
  polynomial_print ( m, o, c, e, title );

  s = - 0.5;
  cout << "\n";
  cout << "  Apply scale factor S = " << s << "\n";
  polynomial_scale ( s, m, o, c, e );

  cout << "\n";
  title = "  S * P(X):";
  polynomial_print ( m, o, c, e, title );

  return;
}
//****************************************************************************80

void polynomial_sort_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_SORT_TEST tests POLYNOMIAL_SORT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c[6] = { 0.0, 9.0, -5.0, - 13.0, 7.0, 11.0 };
  int m = 3;
  int e[6] = { 12, 4, 2, 33, 1, 5 }; 
  int o = 6;
  string title;

  cout << "\n";
  cout << "POLYNOMIAL_SORT_TEST\n";
  cout << "  POLYNOMIAL_SORT sorts a polynomial by exponent index.\n";

  cout << "\n";
  title = "  Unsorted polynomial:";
  polynomial_print ( m, o, c, e, title );

  polynomial_sort ( o, c, e );

  cout << "\n";
  title = "  Sorted polynomial:";
  polynomial_print ( m, o, c, e, title );

  return;
}
//****************************************************************************80

void polynomial_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_VALUE_TEST tests POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c[6] = { 7.0, - 5.0, 9.0, 11.0, 0.0, - 13.0 };
  int m = 3;
  int e[6] = { 1, 2, 4, 5, 12, 33 }; 
  int j;
  int nx = 2;
  int o = 6;
  double *p;
  char title[] = "  P(X) =";
  double x[3*2] = {
     1.0, 2.0, 3.0,
    -2.0, 4.0, 1.0 };

  cout << "\n";
  cout << "POLYNOMIAL_VALUE_TEST\n";
  cout << "  POLYNOMIAL_VALUE evaluates a polynomial.\n";

  cout << "\n";
  polynomial_print ( m, o, c, e, title );

  p = polynomial_value ( m, o, c, e, nx, x );

  cout << "\n";
  for ( j = 0; j < nx; j++ )
  {
    cout << "  P(" << x[0+j*m]
         << "," << x[1+j*m]
         << "," << x[2+j*m]
         << ") = " << p[j] << "\n";
  }

  delete [] p;

  return;
}

