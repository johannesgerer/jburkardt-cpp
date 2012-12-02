# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "divdif.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DIVDIF_PRB.
//
//  Discussion:
//
//    DIVDIF_PRB tests the DIVDIF library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "DIVDIF_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the DIVDIF library.\n";
 
  test01 ( );
  test02 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );
  test11 ( );
  test12 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DIVDIF_PRB\n";
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
//    TEST01 tests DIF_ANTIDERIV, DIF_APPEND, DIF_DERIV_TABLE, DIF_PRINT, DIF_SHIFT_ZERO, DIF_VAL;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 October 2011
//
//  Author:
//
//    John Burkardt
//
{
# define MAXTAB 10

  double diftab[MAXTAB];
  double diftab2[MAXTAB];
  double diftab3[MAXTAB];
  int i;
  int ntab;
  int ntab2;
  int ntab3;
  double xtab[MAXTAB];
  double xtab2[MAXTAB];
  double xtab3[MAXTAB];
  double xval;
  double ytab[MAXTAB];
  double yval;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  DATA_TO_DIF_DISPLAY sets up a difference table\n";
  cout << "  and displays intermediate calculations;\n";
  cout << "  DIF_APPEND appends a new data point;\n";
  cout << "  DIF_ANTIDERIV computes the antiderivative;\n";
  cout << "  DIF_DERIV_TABLE computes the derivative;\n";
  cout << "  DIF_SHIFT_ZERO shifts all the abscissas to 0;\n";
  cout << "  DIF_VAL evaluates at a point;\n";
  cout << "\n";
//
//  Set XTAB, YTAB to X, X^2.
//
  ntab = 4;
  
  for ( i = 0; i < ntab; i++ )
  {
    xtab[i] = ( double ) ( i + 1 );
    ytab[i] = xtab[i] * xtab[i];
  }

  data_to_dif_display ( ntab, xtab, ytab, diftab );

  dif_print ( ntab, xtab, diftab, "  The divided difference polynomial:" );
//
//  Add (5,25) to the table.
//
  cout << "\n";
  cout << "  DIF_APPEND can add the data (5,25) to the table.\n";
  cout << "\n";

  xval = 5.0;
  yval = 25.0;

  dif_append ( ntab, xtab, diftab, xval, yval, &ntab, xtab, diftab );

  dif_print ( ntab, xtab, diftab, 
    "  The updated divided difference polynomial:" );
//
//  Evaluate the polynomial at 2.5.
//
  cout << "\n";
  cout << "  DIF_VAL can evaluate the table at a point.\n";
  cout << "\n";

  xval = 2.5;

  yval = dif_val ( ntab, xtab, diftab, xval );

  cout << "\n";
  cout << "  DIF_VAL reports P(" << xval << ") = " << yval << "\n";
//
//  Shift the base to zero.
//
  dif_shift_zero ( ntab, xtab, diftab );

  dif_print ( ntab, xtab, diftab, 
    "  The divided difference table after DIF_SHIFT_ZERO:" );
//
//  Compute a table for the derivative.
//
  dif_deriv_table ( ntab, xtab, diftab, &ntab2, xtab2, diftab2 );

  dif_print ( ntab2, xtab2, diftab2,
    "  The divided difference table for the derivative:" );

  yval = dif_val ( ntab2, xtab2, diftab2, xval );

  cout << "\n";
  cout << "  DIF_VAL reports P'(" << xval << ") = " << yval << "\n";
//
//  Compute the antiderivative.
//
  dif_antideriv ( ntab, xtab, diftab, &ntab3, xtab3, diftab3 );

  dif_print ( ntab3, xtab3, diftab3,
    "  The divided difference table for the antiderivative:" );

  yval = dif_val ( ntab3, xtab3, diftab3, xval );

  cout << "\n";
  cout << "  DIF_VAL reports (Anti)P(" << xval << ") = " << yval << "\n";

  return;
# undef MAXTAB
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests DATA_TO_DIF and DIF_VAL.
//
//  Discussion:
//
//    This routine demonstrates how divided difference approximation improves with N.
//
//    Evaluate these polynomials at 2.5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXTAB 8

  double diftab[MAXTAB];
  double error;
  int j;
  int ntab;
  double true_value;
  double xtab[MAXTAB];
  double xval;
  double ytab[MAXTAB];
  double yval;
 
  cout << "\n";
  cout << "TEST02\n";
  cout << "  Approximate Y = EXP(X) using orders 1 to " << MAXTAB << ".\n";

  cout << "\n";
  cout << "  Original data:\n";
  cout << "\n";
  cout << "       X          Y\n";
  cout << "\n";
  for ( j = 0; j < MAXTAB; j++ )
  {
    xtab[j] = ( double ) j;
    ytab[j] = exp ( xtab[j] );

    cout                        << "  " 
         << setw(10) << xtab[j] << "  " 
         << setw(10) << ytab[j] << "\n";
  }

  xval = 2.5;
  true_value = exp ( xval );
  cout << "\n";
  cout << "  Evaluate at X = " << xval << " where EXP(X) = " 
    << true_value << "\n";
  cout << "\n";
  cout << "  Order  Approximate Y     Error\n";
  cout << "\n";

  for ( ntab = 1; ntab <= MAXTAB; ntab++ )
  {
 
    for ( j = 0; j < ntab; j++ )
    {
      xtab[j] = ( double ) j;
      ytab[j] = exp ( xtab[j] );
    }

    data_to_dif ( ntab, xtab, ytab, diftab );

    yval = dif_val ( ntab, xtab, diftab, xval );
 
    error = yval - true_value;

    cout                      << "  " 
         << setw(6)  << ntab  << "  " 
         << setw(10) << yval  << "  " 
         << setw(10) << error << "\n";
 
  }
 
  return;
# undef MAXTAB
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests DIF_BASIS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NTAB 5

  double *diftab;
  int i;
  int j;
  int nstep = 9;
  double *pointer;
  int set = 3;
  double xhi;
  double xlo;
  double xtab[NTAB];
  double xval;
  double yval;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  DIF_BASIS computes Lagrange basis polynomials\n";
  cout << "  in difference form.\n";
  cout << "\n";
//
//  Set the base points.
//
  r8vec_indicator ( NTAB, xtab );

  r8vec_print ( NTAB, xtab, "  The base points:" );
//
//  Get the difference tables for the basis polynomials and print them.
//
  diftab = new double[NTAB*NTAB];

  dif_basis ( NTAB, xtab, diftab );

  cout << "\n";
  cout << "  The table of difference vectors defining the basis polynomials.\n";
  cout << "  Each ROW represents a polynomial.\n";
  cout << "\n";

  pointer = diftab;

  for ( i = 0; i < NTAB; i++ )
  {
    cout << "  ";
    for ( j = 0; j < NTAB; j++ )
    {
      cout << setw(10) << *pointer << "  ";
      pointer++;
    }
    cout << "\n";
  }
//
//  Evaluate basis polynomial 3 at a set of points.
//
  cout << "\n";
  cout << "  Evaluate basis polynomial #" << set << " at a set of points.\n";
  cout << "\n";
  cout << "      X        Y\n";
  cout << "\n";

  xhi = ( double ) NTAB;
  xlo = 1.0;
//
//  Advance pointer to beginning of data for basis polynomial SET.
//
  pointer = diftab; 
  for ( i = 1; i < set; i++ )
  {
    for ( j = 1; j <= NTAB; j++ )
    {
      pointer++;
    }
  }

  for ( i = 1; i <= nstep; i++ )
  {
 
    xval = ( ( double ) ( nstep - i     ) * xlo 
           + ( double ) (         i - 1 ) * xhi ) 
           / ( double ) ( nstep     - 1 );
 
    yval = dif_val ( NTAB, xtab, pointer, xval );

    cout                     << "  "
         << setw(10) << xval << "  " 
         << setw(10) << yval << "\n"; 
 
  }
 
  delete [] diftab;

  return;
# undef NTAB
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests POLY_BASIS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NTAB 5

  int i;
  int j;
  int nstep = 9;
  double *pointer;
  double *polcof;
  int set = 3;
  double xhi;
  double xlo;
  double xtab[NTAB];
  double xval;
  double yval;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  POLY_BASIS computes Lagrange basis polynomials\n";
  cout << "  in standard form.\n";
  cout << "\n";

  polcof = new double[NTAB*NTAB];
//
//  Set the base points.
//
  r8vec_indicator ( NTAB, xtab );
//
//  Get the difference tables for the basis polynomials and print them.
//
  r8poly_basis ( NTAB, xtab, polcof );

  cout << "\n";
  cout << "  The table of difference vectors defining the basis polynomials.\n";
  cout << "  Each ROW represents a polynomial.\n";
  cout << "\n";

  pointer = polcof;

  for ( i = 0; i < NTAB; i++ )
  {
    for ( j = 0; j < NTAB; j++ )
    {
      cout << setw(10) << *pointer << "  ";
      pointer++;
    }
    cout << "\n";
  }
//
//  Advance the pointer to the beginning of the data for basis polynomial SET.
//
  pointer = polcof; 
  for ( i = 1; i < set; i++ )
  {
    for ( j = 1; j <= NTAB; j++ )
    {
      pointer++;
    }
  }
//
//  Print basis polynomial SET in polynomial form.
//
  r8poly_print ( NTAB, pointer, "  One basis polynomial in standard form:" );
//
//  Evaluate basis polynoimial SET at a set of points.
//
  cout << "\n";
  cout << "  Evaluate basis polynomial #" << set << " at a set of points.\n";
  cout << "\n";
  cout << "      X        Y\n";
  cout << "\n";

  xhi = ( double ) NTAB;
  xlo = 1.0;

  for ( i = 1; i <= nstep; i++ )
  {
 
    xval = ( ( double ) ( nstep - i     ) * xlo 
           + ( double ) (         i - 1 ) * xhi ) 
           / ( double ) ( nstep     - 1 );
 
    yval = r8poly_val_horner ( NTAB, pointer, xval );
 
    cout                     << "  "
         << setw(10) << xval << "  " 
         << setw(10) << yval << "\n";
 
  }
 
  delete [] polcof;

  return;
# undef NTAB
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests DATA_TO_DIF_DISPLAY, DIF_PRINT, DIF_SHIFT_ZERO, DIF_TO_R8POLY;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXTAB 10

  double c[MAXTAB];
  double diftab1[MAXTAB];
  double diftab2[MAXTAB];
  int i;
  int ntab;
  double x;
  double xtab1[MAXTAB];
  double xtab2[MAXTAB];
  double ytab1[MAXTAB];
  double ytab2[MAXTAB];

  cout << "\n";
  cout << "TEST06\n";
  cout << "  DIF_TO_R8POLY converts a difference table to a polynomial;\n";
  cout << "  DIF_SHIFT_ZERO shifts a divided difference table to 0;\n";
  cout << "\n";
  cout << "  These are equivalent operations\n";
  cout << "\n";
//
//  Set XTAB, YTAB to X, F(X)
//
  ntab = 4;
  for ( i = 0; i < ntab; i++ )
  {
    x = ( double ) ( i + 1 );
    xtab1[i] = x;
    xtab2[i] = x;
    ytab1[i] = - 4.0 + x * ( 3.0 + x * ( - 2.0 + x ) );
    ytab2[i] = ytab1[i];
  }
//
//  Compute and display the finite difference table.
//
  data_to_dif_display ( ntab, xtab1, ytab1, diftab1 );

  data_to_dif_display ( ntab, xtab2, ytab2, diftab2 );
//
//  Examine the corresponding polynomial.
//
  dif_print ( ntab, xtab1, diftab1, "  The divided difference table:" );
//
//  Shift to zero using DIF_SHIFT_ZERO.
//
  dif_shift_zero ( ntab, xtab1, diftab1 );

  r8poly_print ( ntab, diftab1, "  The polynomial using DIF_SHIFT_ZERO:" );
//
//  Shift to zero using DIF_TO_R8POLY.
//
  dif_to_r8poly ( ntab, xtab2, diftab2, c );

  r8poly_print ( ntab, c, "  The polynomial using DIF_TO_R8POLY:" );

  return;

# undef MAXTAB
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests NCC_RULE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NORDER 8

  int i;
  double weight[NORDER];
  double xtab[NORDER];

  cout << "\n";
  cout << "TEST07\n";
  cout << "  NCC_RULE computes closed Newton Cotes formulas;\n";
  cout << "\n";

  ncc_rule ( NORDER, xtab, weight );

  cout << "\n";
  cout << "  Newton-Cotes Closed Quadrature Rule:\n";
  cout << "\n";
  cout << "    Abscissa       Weight\n";
  cout << "\n";

  for ( i = 0; i < NORDER; i++ )
  {
    cout                          << "  "
         << setw(6)  << i+1       << "  " 
         << setw(10) << xtab[i]   << "  "
         << setw(10) << weight[i] << "\n";
  }

  return;
# undef NORDER
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests NCO_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NORDER 8

  int i;
  double weight[NORDER];
  double xtab[NORDER];

  cout << "\n";
  cout << "TEST08\n";
  cout << "  NCO_RULE computes open Newton Cotes formulas.\n";
  cout << "\n";

  nco_rule ( NORDER, xtab, weight );

  cout << "\n";
  cout << "  Newton-Cotes Open Quadrature Rule:\n";
  cout << "\n";
  cout << "    Abscissa       Weight\n";
  cout << "\n";

  for ( i = 0; i < NORDER; i++ )
  {
    cout                          << "  "
         << setw(6)  << i+1       << "  " 
         << setw(10) << xtab[i]   << "  "
         << setw(10) << weight[i] << "\n";
  }

  return;
# undef NORDER
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests R8POLY_ANT_COF, R8POLY_ANT_VAL, R8POLY_DER_COF, R8POLY_DER_VAL, R8POLY_PRINT, R8POLY_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  double poly_cof[N];
  double poly_cof2[N+1];
  double poly_cof3[N-1];
  double xval;
  double yval;
  double yval2;
  double yval3;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  R8POLY_ANT_COF computes the coefficients of the\n";
  cout << "  antiderivative of a polynomial;\n";
  cout << "  R8POLY_ANT_VAL evaluates the antiderivative of\n";
  cout << "  a polynomial;\n";
  cout << "  R8POLY_DER_COF computes the coefficients of the\n";
  cout << "  derivative of a polynomial;\n";
  cout << "  R8POLY_DER_VAL evaluates the derivative of\n";
  cout << "  a polynomial;\n";
  cout << "  R8POLY_PRINT prints a polynomial;\n";
  cout << "  R8POLY_VAL evaluates a polynomial.\n";

  for ( i = 0; i <= N; i++ )
  {
    poly_cof[i] = ( double ) ( i + 1 );
  }

  r8poly_print ( N, poly_cof, "  The initial polynomial:" );

  r8poly_ant_cof ( N, poly_cof, poly_cof2 );

  r8poly_print ( N+1, poly_cof2, "  The antiderivative polynomial:" );

  r8poly_der_cof ( N, poly_cof, poly_cof3 );

  r8poly_print ( N-1, poly_cof3, "  The derivative polynomial:" );

  cout << "\n";
  cout << "  Evaluate the polynomial, antiderivative and\n";
  cout << "  derivative, using only the original polynomial\n";
  cout << "  coefficients:\n";
  cout << "\n";
  cout << "  X   P(X)   Anti_P(X)     P'(X)\n";
  cout << "\n";

  for ( i = 0; i <= 2; i++ )
  {

    xval = ( double ) i;


    yval = r8poly_val_horner ( N, poly_cof, xval );

    yval2 = r8poly_ant_val ( N, poly_cof, xval );

    yval3 = r8poly_der_val ( N, poly_cof, xval );

    cout                      << "  " 
         << setw(10) << xval  << "  " 
         << setw(10) << yval  << "  " 
         << setw(10) << yval2 << "  " 
         << setw(10) << yval3 << "\n";

  }

  return;
# undef N
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests ROOTS_TO_DIF and DIF_TO_R8POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXROOTS 4

  double c[MAXROOTS+1];
  double diftab[MAXROOTS+1];
  int nroots;
  int ntab;
  double roots[MAXROOTS];
  double xtab[MAXROOTS+1];

  cout << "\n";
  cout << "TEST10\n";
  cout << "  ROOTS_TO_DIF computes a divided difference\n";
  cout << "  polynomial with the given roots;\n";
  cout << "  DIF_TO_R8POLY converts it to a standard form polynomial.\n";
  cout << "\n";

  nroots = 1;
  roots[0] = 3.0;
  r8vec_print ( nroots, roots, "  The roots:" );

  roots_to_dif ( nroots, roots, &ntab, xtab, diftab );
  dif_to_r8poly ( ntab, xtab, diftab, c );
  r8poly_print ( ntab, c, "  The polynomial:" );
 
  nroots = 2;
  roots[0] = 3.0;
  roots[1] = 1.0;
  r8vec_print ( nroots, roots, "  The roots:" );

  roots_to_dif ( nroots, roots, &ntab, xtab, diftab );
  dif_to_r8poly ( ntab, xtab, diftab, c );
  r8poly_print ( ntab, c, "  The polynomial:" );

  nroots = 3;
  roots[0] = 3.0;
  roots[1] = 1.0;
  roots[2] = 2.0;
  r8vec_print ( nroots, roots, "  The roots:" );

  roots_to_dif ( nroots, roots, &ntab, xtab, diftab );
  dif_to_r8poly ( ntab, xtab, diftab, c );
  r8poly_print ( ntab, c, "  The polynomial:" );

  nroots = 4;
  roots[0] = 3.0;
  roots[1] = 1.0;
  roots[2] = 2.0;
  roots[3] = 4.0;
  r8vec_print ( nroots, roots, "  The roots:" );

  roots_to_dif ( nroots, roots, &ntab, xtab, diftab );
  dif_to_r8poly ( ntab, xtab, diftab, c );
  r8poly_print ( ntab, c, "  The polynomial:" );

  return;
# undef MAXTAB
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests ROOTS_TO_R8POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXROOT 5

  double c[MAXROOT+1];
  int nc;
  int nroot;
  double roots[MAXROOT];

  cout << "\n";
  cout << "TEST11\n";
  cout << "  ROOTS_TO_R8POLY computes polynomial coefficients from roots.\n";
  cout << "\n";
 
  nroot = 1;
  roots[0] = 3.0;
  r8vec_print ( nroot, roots, "  The roots:" );

  roots_to_r8poly ( nroot, roots, &nc, c );

  r8poly_print ( nc, c, "  The polynomial:" );
 
  nroot = 2;
  roots[0] = 3.0;
  roots[1] = 1.0;
  r8vec_print ( nroot, roots, "  The roots:" );

  roots_to_r8poly ( nroot, roots, &nc, c );

  r8poly_print ( nc, c, "  The polynomial:" );

  nroot = 3;
  roots[0] = 3.0;
  roots[1] = 1.0;
  roots[2] = 2.0;
  r8vec_print ( nroot, roots, "  The roots:" );

  roots_to_r8poly ( nroot, roots, &nc, c );

  r8poly_print ( nc, c, "  The polynomial:" );

  nroot = 4;
  roots[0] = 3.0;
  roots[1] = 1.0;
  roots[2] = 2.0;
  roots[3] = 4.0;
  r8vec_print ( nroot, roots, "  The roots:" );

  roots_to_r8poly ( nroot, roots, &nc, c );

  r8poly_print ( nc, c, "  The polynomial:" );

  return;

# undef MAXROOT
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests R8POLY_SHIFT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int i;
  double poly_cof[N];
  double scale;
  double shift;

  scale = 2.0;
  shift = +3.0;
  poly_cof[0] = +6.0;
  poly_cof[1] = -1.0;
  poly_cof[2] =  2.0;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  R8POLY_SHIFT shifts polynomial coefficients.\n";
  cout << "\n";
  cout << "  Polynomial coefficients for argument X\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout                            << "  "
         << setw(6)  << i           << "  " 
         << setw(10) << poly_cof[i] << "\n";
  }

  r8poly_shift ( scale, shift, N, poly_cof );

  cout << "\n";
  cout << "  SCALE = " << scale << "\n";
  cout << "  SHIFT = " << shift << "\n";
  cout << "\n";
  cout << "  Polynomial coefficients for argument Z = SCALE * X + SHIFT:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout                            << "  "
         << setw(6)  << i           << "  " 
         << setw(10) << poly_cof[i] << "\n";
  }

  return;

# undef N
}
