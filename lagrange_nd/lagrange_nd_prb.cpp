# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <stdio.h>

using namespace std;

# include "lagrange_nd.hpp"

int main ( );

void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( int option );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LAGRANGE_ND_PRB.
//
//  Discussion:
//
//    LAGRANGE_ND_PRB tests the LAGRANGE_ND library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int option;

  timestamp ( );
  cout << "\n";
  cout << "LAGRANGE_ND_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LAGRANGE_ND library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );

  option = 0;
  test11 ( option );

  option = 1;
  test11 ( option );
//
//  Terminate.
//
  cout << "\n";
  cout << "LAGRANGE_ND_PRB\n";
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
//    TEST01 tests MONO_BETWEEN_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int n1;
  int n2;
  int v;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  MONO_BETWEEN_ENUM can enumerate the number of monomials\n";
  cout << "  in D variables, of total degree between N1 and N2.\n";

  d = 3;
  cout << "\n";
  cout << "  Using spatial dimension D = " << d << "\n";
  cout << "\n";
  cout << "   N2:";
  for ( n2 = 0; n2 <= 8; n2++ )
  {
    cout << "  " << setw(4) << n2;
  }
  cout << "\n";
  cout << "  N1 +------------------------------------------------------\n";
  for ( n1 = 0; n1 <= 8; n1++ )
  {
    cout << "  " << setw(2) << n1 << " |";
    for ( n2 = 0; n2 <= 8; n2++ )
    {
      v = mono_between_enum ( d, n1, n2 );
      cout << "  " << setw(4) << v;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests MONO_TOTAL_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int n;
  int v;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  MONO_TOTAL_ENUM can enumerate the number of monomials\n";
  cout << "  in D variables, of total degree N.\n";

  cout << "\n";
  cout << "    N:";
  for ( n = 0; n <= 8; n++ )
  {
    cout << setw(4) << n;
  }
  cout << "\n";
  cout << "   D +------------------------------------------------------\n";
  for ( d = 1; d <= 8; d++ )
  {
    cout << "  " << setw(2) << d << " |";
    for ( n = 0; n <= 8; n++ )
    {
      v = mono_total_enum ( d, n );
      cout << "  " << setw(4) << v;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests MONO_UPTO_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int n;
  int v;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  MONO_UPTO_ENUM can enumerate the number of monomials\n";
  cout << "  in D variables, of total degree 0 up to N.\n";

  cout << "\n";
  cout << "    N:";
  for ( n = 0; n <= 8; n++ )
  {
    cout << "  " << setw(4) << n;
  }
  cout << "\n";
  cout << "   D +------------------------------------------------------\n";
  for ( d = 1; d <= 8; d++ )
  {
    cout << "  " << setw(2) << d << " |";
    for ( n = 0; n <= 8; n++ )
    {
      v = mono_upto_enum ( d, n );
      cout << " " << setw(5) << v;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests MONO_BETWEEN_NEXT_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int d = 3;
  int i;
  int j;
  int n1;
  int n2;
  int x[3];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  MONO_BETWEEN_NEXT_GRLEX can list the monomials\n";
  cout << "  in D variables, of total degree N between N1 and N2,\n";
  cout << "  one at a time.\n";
  cout << "\n";
  cout << "  We start the process with (0,0,...,0,N1).\n";
  cout << "  The process ends with (N2,0,...,0,0)\n";

  n1 = 2;
  n2 = 3;

  cout << "\n";
  cout << "  Let D =  " << d << "\n";
  cout << "      N1 = " << n1 << "\n";
  cout << "      N2 = " << n2 << "\n";
  cout << "\n";

  x[0] = 0;
  x[1] = 0;
  x[2] = n1;

  j = 1;

  for ( ; ; )
  {
    cout << "  " << setw(2) << j << ":";
    for ( i = 0; i < d; i++ )
    {
      cout << "  " << setw(1) << x[i];
    }
    cout << "\n";

    if ( x[0] == n2 )
    {
      break;
    }

    mono_between_next_grlex ( d, n1, n2, x );
    j = j + 1;
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests LAGRANGE_ND_COMPLETE in 1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int d;
  int *e;
  double error;
  int i;
  int j;
  char label[80];
  int n;
  int nd;
  int o;
  double *pc;
  int *pe;
  int *po;
  int r;
  double *v;
  double *value;
  double xd[5] = { 0.0, 1.0, 2.0, 3.0, 4.0 };

  cout << "\n";
  cout << "TEST05\n";
  cout << "  LAGRANGE_COMPLETE determines\n";
  cout << "  the Lagrange interpolating polynomials L(x)\n";
  cout << "  for ND points in D dimensions, assuming that\n";
  cout << "  the number of points exactly coincides with\n";
  cout << "  R = Pi(D,N), the number of monomials of degree N or less\n";
  cout << "\n";
  cout << "  As a special demonstration, this code runs in 1D\n";

  nd = 5;

  d = 1;
  n = 4;
  r = mono_upto_enum ( d, n );

  pc = new double[nd*r];
  pe = new int[nd*r];
  po = new int[nd];

  c = new double[r];
  e = new int[r];

  cout << "\n";
  cout << "  Spatial dimension D = " << d << "\n";
  cout << "  Maximum degree N = " << n << "\n";
  cout << "  Number of monomials R = " << r << "\n";
  cout << "  Number of data points ND = " << nd << "\n";

  r8mat_transpose_print ( d, nd, xd, "  Data points XD:" );

  lagrange_complete ( d, n, r, nd, xd, po, pc, pe );
//
//  Print the polynomials.
//
  cout << "\n";
  cout << "  Lagrange polynomials for XD data points:\n";
  cout << "\n";

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    sprintf ( label, "  P(%d)(x) =", i );
    polynomial_print ( d, o, c, e, label );
  }
//
//  Evaluate the polynomials at XD.
//
  value = new double[nd*nd];

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    v = polynomial_value ( d, o, c, e, nd, xd );

    for ( j = 0; j < nd; j++ )
    {
      value[i+j*nd] = v[j];
    }
    delete [] v;
  }

  error = r8mat_is_identity ( nd, value );
  cout << "\n";
  cout << "  Frobenius norm of Lagrange matrix error = " << error << "\n";
//
//  Free memory.
//
  delete [] c;
  delete [] e;
  delete [] pc;
  delete [] pe;
  delete [] po;
  delete [] value;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests LAGRANGE_ND_COMPLETE in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int d;
  int *e;
  double error;
  int i;
  int j;
  char label[80];
  int n;
  int nd;
  int o;
  double *pc;
  int *pe;
  int *po;
  int r;
  double *v;
  double *value;
  double xd[2*6] = {
    0.0,  0.0, 
    1.0,  0.0, 
    2.0,  0.0, 
    0.0,  1.0, 
    1.0,  1.0, 
    0.0,  2.0 };

  cout << "\n";
  cout << "TEST06\n";
  cout << "  LAGRANGE_COMPLETE determines\n";
  cout << "  the Lagrange interpolating polynomials L(x)\n";
  cout << "  for ND points in D dimensions, assuming that\n";
  cout << "  the number of points exactly coincides with\n";
  cout << "  R = Pi(D,N), the number of monomials of degree N or less\n";
  cout << "\n";
  cout << "  The data points are the grid nodes of a triangle.\n";

  nd = 6;

  d = 2;
  n = 2;
  r = mono_upto_enum ( d, n );

  pc = new double[nd*r];
  pe = new int[nd*r];
  po = new int[nd];

  c = new double[r];
  e = new int[r];

  cout << "\n";
  cout << "  Spatial dimension D = " << d << "\n";
  cout << "  Maximum degree N = " << n << "\n";
  cout << "  Number of monomials R = " << r << "\n";
  cout << "  Number of data points ND = " << nd << "\n";

  r8mat_transpose_print ( d, nd, xd, "  Data points XD:" );

  lagrange_complete ( d, n, r, nd, xd, po, pc, pe );
//
//  Print the polynomials.
//
  cout << "\n";
  cout << "  Lagrange polynomials for XD data points:\n";
  cout << "\n";

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    sprintf ( label, "  P(%d)(x) =", i );
    polynomial_print ( d, o, c, e, label );
  }
//
//  Evaluate the polynomials at XD.
//
  value = new double[nd*nd];

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    v = polynomial_value ( d, o, c, e, nd, xd );

    for ( j = 0; j < nd; j++ )
    {
      value[i+j*nd] = v[j];
    }
    delete [] v;
  }

  error = r8mat_is_identity ( nd, value );
  cout << "\n";
  cout << "  Frobenius norm of Lagrange matrix error = " << error << "\n";
//
//  Free memory.
//
  delete [] c;
  delete [] e;
  delete [] pc;
  delete [] pe;
  delete [] po;
  delete [] value;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests LAGRANGE_ND_COMPLETE in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int d;
  int *e;
  double error;
  int i;
  int j;
  char label[80];
  int n;
  int nd;
  int o;
  double *pc;
  int *pe;
  int *po;
  int r;
  double *v;
  double *value;
  double xd[3*10] = {
    0.0,  0.0,  0.0, 
    1.0,  0.0,  0.0, 
    2.0,  0.0,  0.0, 
    0.0,  1.0,  0.0, 
    1.0,  1.0,  0.0, 
    0.0,  2.0,  0.0, 
    0.0,  0.0,  1.0, 
    1.0,  0.0,  1.0, 
    0.0,  1.0,  1.0, 
    0.0,  0.0,  2.0 };

  cout << "\n";
  cout << "TEST07\n";
  cout << "  LAGRANGE_COMPLETE determines\n";
  cout << "  the Lagrange interpolating polynomials L(x)\n";
  cout << "  for ND points in D dimensions, assuming that\n";
  cout << "  the number of points exactly coincides with\n";
  cout << "  R = Pi(D,N), the number of monomials of degree N or less\n";
  cout << "\n";
  cout << "  The data points are the grid nodes of a tetrahedron.\n";

  nd = 10;

  d = 3;
  n = 2;
  r = mono_upto_enum ( d, n );

  pc = new double[nd*r];
  pe = new int[nd*r];
  po = new int[nd];

  c = new double[r];
  e = new int[r];

  cout << "\n";
  cout << "  Spatial dimension D = " << d << "\n";
  cout << "  Maximum degree N = " << n << "\n";
  cout << "  Number of monomials R = " << r << "\n";
  cout << "  Number of data points ND = " << nd << "\n";

  r8mat_transpose_print ( d, nd, xd, "  Data points XD:" );

  lagrange_complete ( d, n, r, nd, xd, po, pc, pe );
//
//  Print the polynomials.
//
  cout << "\n";
  cout << "  Lagrange polynomials for XD data points:\n";
  cout << "\n";

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    sprintf ( label, "  P(%d)(x) =", i );
    polynomial_print ( d, o, c, e, label );
  }
//
//  Evaluate the polynomials at XD.
//
  value = new double[nd*nd];

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    v = polynomial_value ( d, o, c, e, nd, xd );

    for ( j = 0; j < nd; j++ )
    {
      value[i+j*nd] = v[j];
    }
    delete [] v;
  }

  error = r8mat_is_identity ( nd, value );
  cout << "\n";
  cout << "  Frobenius norm of Lagrange matrix error = " << error << "\n";
//
//  Free memory.
//
  delete [] c;
  delete [] e;
  delete [] pc;
  delete [] pe;
  delete [] po;
  delete [] value;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests LAGRANGE_PARTIAL in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int d;
  int *e;
  double error;
  int i;
  int j;
  char label[80];
  int n;
  int nd;
  int o;
  double *pc;
  int *pe;
  int *po;
  int r;
  double *v;
  double *value;
  double xd[2*5] = {
     0.0,  0.0, 
    -1.0,  0.0, 
     1.0,  0.0, 
     0.0, -1.0, 
     0.0,  1.0 };

  cout << "\n";
  cout << "TEST08\n";
  cout << "  LAGRANGE_PARTIAL determines\n";
  cout << "  the Lagrange interpolating polynomials L(x)\n";
  cout << "  for ND points in D dimensions, assuming that\n";
  cout << "  the number of points is less than or equal to\n";
  cout << "  R = Pi(D,N), the number of monomials of degree N or less\n";
  cout << "\n";
  cout << "  For this example, the data points are the same as those\n";
  cout << "  used by the level 1 Clenshaw Curtis sparse grid in 2D.\n";

  nd = 5;
 
  d = 2;
  n = 2;
  r = mono_upto_enum ( d, n );

  pc = new double[nd*r];
  pe = new int[nd*r];
  po = new int[nd];

  c = new double[r];
  e = new int[r];

  cout << "\n";
  cout << "  Spatial dimension D = " << d << "\n";
  cout << "  Maximum degree N = " << n << "\n";
  cout << "  Number of monomials R = " << r << "\n";
  cout << "  Number of data points ND = " << nd << "\n";

  r8mat_transpose_print ( d, nd, xd, "  Data points XD:" );

  lagrange_partial ( d, n, r, nd, xd, po, pc, pe );
//
//
//  Print the polynomials.
//
  cout << "\n";
  cout << "  Lagrange polynomials for XD data points:\n";
  cout << "\n";

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    sprintf ( label, "  P(%d)(x) =", i );
    polynomial_print ( d, o, c, e, label );
  }
//
//  Evaluate the polynomials at XD.
//
  value = new double[nd*nd];

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    v = polynomial_value ( d, o, c, e, nd, xd );

    for ( j = 0; j < nd; j++ )
    {
      value[i+j*nd] = v[j];
    }
    delete [] v;
  }

  error = r8mat_is_identity ( nd, value );
  cout << "\n";
  cout << "  Frobenius norm of Lagrange matrix error = " << error << "\n";
//
//  Free memory.
//
  delete [] c;
  delete [] e;
  delete [] pc;
  delete [] pe;
  delete [] po;
  delete [] value;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests LAGRANGE_PARTIAL in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int d;
  int *e;
  double error;
  int i;
  int j;
  char label[80];
  int n;
  int nd;
  int o;
  double *pc;
  int *pe;
  int *po;
  int r;
  const double sq2h = sqrt ( 2.0 ) / 2.0;
  double *v;
  double *value;
  double xd[3*25] = {
      0.0,   0.0,   0.0, 
     -1.0,   0.0,   0.0, 
      1.0,   0.0,   0.0, 
      0.0,  -1.0,   0.0, 
      0.0,   1.0,   0.0, 
      0.0,   0.0,  -1.0, 
      0.0,   0.0,   1.0, 
     -0.707106781187,  0.0,   0.0, 
      0.707106781187,  0.0,   0.0, 
     -1.0,  -1.0,   0.0, 
      1.0,  -1.0,   0.0, 
     -1.0,   1.0,   0.0, 
      1.0,   1.0,   0.0, 
      0.0,  -0.707106781187,  0.0, 
      0.0,   0.707106781187,  0.0, 
     -1.0,   0.0,  -1.0, 
      1.0,   0.0,  -1.0, 
     -1.0,   0.0,   1.0, 
      1.0,   0.0,   1.0, 
      0.0,  -1.0,  -1.0, 
      0.0,   1.0,  -1.0, 
      0.0,  -1.0,   1.0, 
      0.0,   1.0,   1.0, 
      0.0,   0.0,  -0.707106781187, 
      0.0,   0.0,   0.707106781187 };

  cout << "\n";
  cout << "TEST09\n";
  cout << "  LAGRANGE_PARTIAL determines\n";
  cout << "  the Lagrange interpolating polynomials L(x)\n";
  cout << "  for ND points in D dimensions, assuming that\n";
  cout << "  the number of points is less than or equal to\n";
  cout << "  R = Pi(D,N), the number of monomials of degree N or less\n";
  cout << "\n";
  cout << "  For this example, the data points are the same as those\n";
  cout << "  used by the level 2 Clenshaw Curtis sparse grid in 3D.\n";

  nd = 25;

  d = 3;
  n = 4;
  r = mono_upto_enum ( d, n );

  pc = new double[nd*r];
  pe = new int[nd*r];
  po = new int[nd];

  c = new double[r];
  e = new int[r];

  cout << "\n";
  cout << "  Spatial dimension D = " << d << "\n";
  cout << "  Maximum degree N = " << n << "\n";
  cout << "  Number of monomials R = " << r << "\n";
  cout << "  Number of data points ND = " << nd << "\n";

  r8mat_transpose_print ( d, nd, xd, "  Data points XD:" );

  lagrange_partial ( d, n, r, nd, xd, po, pc, pe );
//
//
//  Print the polynomials.
//
  cout << "\n";
  cout << "  Lagrange polynomials for XD data points:\n";
  cout << "\n";

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    sprintf ( label, "  P(%d)(x) =", i );
    polynomial_print ( d, o, c, e, label );
  }
//
//  Evaluate the polynomials at XD.
//
  value = new double[nd*nd];

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    v = polynomial_value ( d, o, c, e, nd, xd );

    for ( j = 0; j < nd; j++ )
    {
      value[i+j*nd] = v[j];
    }
    delete [] v;
  }

  error = r8mat_is_identity ( nd, value );
  cout << "\n";
  cout << "  Frobenius norm of Lagrange matrix error = " << error << "\n";
//
//  Free memory.
//
  delete [] c;
  delete [] e;
  delete [] pc;
  delete [] pe;
  delete [] po;
  delete [] value;

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests LAGRANGE_PARTIAL2 in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int d;
  int *e;
  double error;
  double f;
  int i;
  int j;
  int k;
  char label[80];
  int n;
  int nd;
  int ni;
  int o;
  double *pc;
  double *pd;
  int *pe;
  int pn;
  int *po;
  int r;
  double *v;
  double *value;
  double xd[2*13] = {
    0.0,  0.0, 
   -1.0,  0.0, 
    1.0,  0.0, 
    0.0, -1.0, 
    0.0,  1.0, 
   -1.0,  1.0, 
    1.0,  1.0, 
   -1.0, -1.0,
    1.0, -1.0, 
   -0.5,  0.0, 
    0.0, -0.5, 
    0.0, +0.5, 
    0.5,  0.0 };
  double *xyi;
  double *zi;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  LAGRANGE_PARTIAL2 determines\n";
  cout << "  the Lagrange interpolating polynomials L(x)\n";
  cout << "  for ND points in D dimensions, assuming that\n";
  cout << "  the number of points is less than or equal to\n";
  cout << "  R = Pi(D,N), the number of monomials of degree N or less\n";
  cout << "\n";
  cout << "  For this example, the data points are the same as those\n";
  cout << "  used by the level 2 Clenshaw Curtis sparse grid in 2D.\n";

  nd = 13;
  ni = 11 * 11;
 
  d = 2;
  n = 4;
  r = mono_upto_enum ( d, n );

  pc = new double[nd*r];
  pe = new int[nd*r];
  po = new int[nd];

  c = new double[r];
  e = new int[r];

  cout << "\n";
  cout << "  Spatial dimension D = " << d << "\n";
  cout << "  Maximum degree N = " << n << "\n";
  cout << "  Number of monomials R = " << r << "\n";
  cout << "  Number of data points ND = " << nd << "\n";

  r8mat_transpose_print ( d, nd, xd, "  Data points XD:" );

  lagrange_partial2 ( d, n, r, nd, xd, po, pc, pe );
//
//  Print the polynomials.
//
  cout << "\n";
  cout << "  Lagrange polynomials for XD data points:\n";
  cout << "\n";

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    sprintf ( label, "  P(%d)(x) =", i );
    polynomial_print ( d, o, c, e, label );
  }
//
//  Evaluate the polynomials at XD.
//
  value = new double[nd*nd];

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    v = polynomial_value ( d, o, c, e, nd, xd );

    for ( j = 0; j < nd; j++ )
    {
      value[i+j*nd] = v[j];
    }
    delete [] v;
  }

  error = r8mat_is_identity ( nd, value );
  cout << "\n";
  cout << "  Frobenius norm of Lagrange matrix error = " << error << "\n";
//
//  Evaluate a function at the data points.
//
  pd = new double[nd];
  for ( i = 0; i < nd; i++ )
  {
    pd[i] = sin ( xd[0+i*2] ) * cos ( xd[1+i*2] );
  }
//
//  Compare exact function and interpolant at a grid of points.
//
  xyi = new double[2*ni];

  k = 0;
  for ( j = 1; j <= 11; j++ )
  {
    for ( i = 1; i <= 11; i++ )
    {
      xyi[0+k*2] = ( ( double ) ( 11 - i     ) * ( - 1.0 )   
                   + ( double ) (      i - 1 ) * ( + 1.0 ) ) 
                   / ( double ) ( 11     - 1 );
      xyi[1+k*2] = ( ( double ) ( 11 - j     ) * ( - 1.0 )   
                   + ( double ) (      j - 1 ) * ( + 1.0 ) ) 
                   / ( double ) ( 11     - 1 );
      k = k + 1;
    }
  }

  pn = nd;
  zi = interpolant_value ( d, r, pn, po, pc, pe, pd, ni, xyi );

  error = 0.0;
  for ( k = 0; k < ni; k++ )
  {
    f = sin ( xyi[0+k*2] ) * cos ( xyi[1+k*2] );
    if ( error < fabs ( zi[k] - f ) )
    {
      error = fabs ( zi[k] - f );
    }
  }
  cout << "\n";
  cout << "  Maximum absolute interpolant error on 11x11 grid = " << error << "\n";
//
//  Free memory.
//
  delete [] c;
  delete [] e;
  delete [] pc;
  delete [] pd;
  delete [] pe;
  delete [] po;
  delete [] value;
  delete [] xyi;
  delete [] zi;

  return;
}
//****************************************************************************80

void test11 ( int option )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_ND_TEST11 tests LAGRANGE_PARTIAL3 in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int OPTION, determines the initial basis:
//    0, use monomials, 1, x, y, x^2, xy, y^2, x^3, ...
//    1, use Legendre products, 1, y, x, (3y^2-1)/2, xy, (3^x^2-1), (5y^3-3)/2,...
//
{
  double *c;
  int d;
  int *e;
  double error;
  double f;
  int i;
  int j;
  int k;
  char label[80];
  int n;
  int n2;
  int nd;
  int ni;
  int o;
  double *pc;
  double *pd;
  int *pe;
  int pn;
  int *po;
  int r;
  double *v;
  double *value;
  double xd[2*65] = {
      0.0000000000000000,      0.0000000000000000,
     -1.0000000000000000,      0.0000000000000000,
      1.0000000000000000,      0.0000000000000000,
      0.0000000000000000,     -1.0000000000000000,
      0.0000000000000000,      1.0000000000000000,
     -0.7071067811865475,      0.0000000000000000,
      0.7071067811865476,      0.0000000000000000,
     -1.0000000000000000,     -1.0000000000000000,
      1.0000000000000000,     -1.0000000000000000,
     -1.0000000000000000,      1.0000000000000000,
      1.0000000000000000,      1.0000000000000000,
      0.0000000000000000,     -0.7071067811865475,
      0.0000000000000000,      0.7071067811865476,
     -0.9238795325112867,      0.0000000000000000,
     -0.3826834323650897,      0.0000000000000000,
      0.3826834323650898,      0.0000000000000000,
      0.9238795325112867,      0.0000000000000000,
     -0.7071067811865475,     -1.0000000000000000,
      0.7071067811865476,     -1.0000000000000000,
     -0.7071067811865475,      1.0000000000000000,
      0.7071067811865476,      1.0000000000000000,
     -1.0000000000000000,     -0.7071067811865475,
      1.0000000000000000,     -0.7071067811865475,
     -1.0000000000000000,      0.7071067811865476,
      1.0000000000000000,      0.7071067811865476,
      0.0000000000000000,     -0.9238795325112867,
      0.0000000000000000,     -0.3826834323650897,
      0.0000000000000000,      0.3826834323650898,
      0.0000000000000000,      0.9238795325112867,
     -0.9807852804032304,      0.0000000000000000,
     -0.8314696123025453,      0.0000000000000000,
     -0.5555702330196020,      0.0000000000000000,
     -0.1950903220161282,      0.0000000000000000,
      0.1950903220161283,      0.0000000000000000,
      0.5555702330196023,      0.0000000000000000,
      0.8314696123025452,      0.0000000000000000,
      0.9807852804032304,      0.0000000000000000,
     -0.9238795325112867,     -1.0000000000000000,
     -0.3826834323650897,     -1.0000000000000000,
      0.3826834323650898,     -1.0000000000000000,
      0.9238795325112867,     -1.0000000000000000,
     -0.9238795325112867,      1.0000000000000000,
     -0.3826834323650897,      1.0000000000000000,
      0.3826834323650898,      1.0000000000000000,
      0.9238795325112867,      1.0000000000000000,
     -0.7071067811865475,     -0.7071067811865475,
      0.7071067811865476,     -0.7071067811865475,
     -0.7071067811865475,      0.7071067811865476,
      0.7071067811865476,      0.7071067811865476,
     -1.0000000000000000,     -0.9238795325112867,
      1.0000000000000000,     -0.9238795325112867,
     -1.0000000000000000,     -0.3826834323650897,
      1.0000000000000000,     -0.3826834323650897,
     -1.0000000000000000,      0.3826834323650898,
      1.0000000000000000,      0.3826834323650898,
     -1.0000000000000000,      0.9238795325112867,
      1.0000000000000000,      0.9238795325112867,
      0.0000000000000000,     -0.9807852804032304,
      0.0000000000000000,     -0.8314696123025453,
      0.0000000000000000,     -0.5555702330196020,
      0.0000000000000000,     -0.1950903220161282,
      0.0000000000000000,      0.1950903220161283,
      0.0000000000000000,      0.5555702330196023,
      0.0000000000000000,      0.8314696123025452,
      0.0000000000000000,      0.9807852804032304 };
  double *xyi;
  double *zi;

  cout << "\n";
  cout << "LAGRANGE_ND_TEST11\n";
  cout << "  LAGRANGE_PARTIAL3 determines\n";
  cout << "  the Lagrange interpolating polynomials L(x)\n";
  cout << "  for ND points in D dimensions, assuming that\n";
  cout << "  the number of points is less than or equal to\n";
  cout << "  R = Pi(D,N), the number of monomials of degree N or less\n";
  cout << "\n";
  cout << "  If LAGRANGE_PARTIAL3 determines that the problem is not\n";
  cout << "  well-posed for the given value of N, it increases N\n";
  cout << "  until a suitable value is found.\n";
  cout << "\n";
  cout << "  For this example, the data points are the same as those\n";
  cout << "  used by the level 2 Clenshaw Curtis sparse grid in 2D.\n";

  nd = 65;
  ni = 11 * 11;

  d = 2;
  n = 10;
  po = new int[nd];

  cout << "\n";
  cout << "  Spatial dimension D = " << d << "\n";
  cout << "  Maximum degree N = " << n << "\n";
  cout << "  Number of data points ND = " << nd << "\n";
  cout << "  Monomial/Legendre option OPTION = " << option << "\n";

  r8mat_transpose_print ( d, nd, xd, "  Data points XD:" );

  lagrange_partial3 ( d, n, nd, xd, option, po, &pc, &pe, n2 );

  if ( n < n2 )
  {
    cout << "\n";
    cout << "  LAGRANGE_PARTIAL3 increased N to " << n2 << "\n";
  }

  r = mono_upto_enum ( d, n2 );
  cout << "  Number of monomials R = " << r << "\n";
  c = new double[r];
  e = new int[r];
//
//  Print the polynomials.
//
  cout << "\n";
  cout << "  (First 2) Lagrange polynomials for XD data points:\n";
  cout << "\n";

//for ( i = 0; i < nd; i++ )
  for ( i = 0; i < 2; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    sprintf ( label, "  P(%d)(x) =", i );
    polynomial_print ( d, o, c, e, label );
  }
//
//  Evaluate the polynomials at XD.
//
  value = new double[nd*nd];

  for ( i = 0; i < nd; i++ )
  {
    o = po[i];
    for ( j = 0; j < o; j++ )
    {
      c[j] = pc[i+j*nd];
      e[j] = pe[i+j*nd];
    }
    v = polynomial_value ( d, o, c, e, nd, xd );

    for ( j = 0; j < nd; j++ )
    {
      value[i+j*nd] = v[j];
    }
    delete [] v;
  }

  error = r8mat_is_identity ( nd, value );
  cout << "\n";
  cout << "  Frobenius norm of Lagrange matrix error = " << error << "\n";
//
//  Evaluate a function at the data points.
//
  pd = new double[nd];
  for ( i = 0; i < nd; i++ )
  {
    pd[i] = exp ( xd[0+i*2] * xd[1+i*2] );
  }
//
//  Compare exact function and interpolant at a grid of points.
//
  xyi = new double[2*ni];

  k = 0;
  for ( j = 1; j <= 11; j++ )
  {
    for ( i = 1; i <= 11; i++ )
    {
      xyi[0+k*2] = ( ( double ) ( 11 - i     ) * ( - 1.0 )   
                   + ( double ) (      i - 1 ) * ( + 1.0 ) ) 
                   / ( double ) ( 11     - 1 );
      xyi[1+k*2] = ( ( double ) ( 11 - j     ) * ( - 1.0 )   
                   + ( double ) (      j - 1 ) * ( + 1.0 ) ) 
                   / ( double ) ( 11     - 1 );
      k = k + 1;
    }
  }

  pn = nd;
  zi = interpolant_value ( d, r, pn, po, pc, pe, pd, ni, xyi );

  error = 0.0;
  for ( k = 0; k < ni; k++ )
  {
    f = exp ( xyi[0+k*2] * xyi[1+k*2] );
    if ( error < fabs ( zi[k] - f ) )
    {
      error = fabs ( zi[k] - f );
    }
  }
  cout << "\n";
  cout << "  Maximum absolute interpolant error on 11x11 grid = " << error << "\n";
//
//  Free memory.
//
  delete [] c;
  delete [] e;
  delete [] pc;
  delete [] pd;
  delete [] pe;
  delete [] po;
  delete [] value;
  delete [] xyi;
  delete [] zi;

  return;
}
