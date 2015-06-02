# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "divdif.hpp"

//****************************************************************************80

double *cheby_t_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T_ZERO returns zeroes of the Chebyshev polynomial T(N)(X).
//
//  Discussion:
//
//    The I-th zero of T(N)(X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Output, double CHEBY_T_ZERO[N], the zeroes of T(N)(X).
//
{
  double angle;
  int i;
  static double pi = 3.141592653589793;
  double *z;

  z = new double[n];

  for ( i = 0; i < n; i++ )
  {
    angle = ( double ) ( 2 * i + 1 ) * pi / ( double ) ( 2 * n );
    z[i] = cos ( angle );
  }
  return z;
}
//****************************************************************************80

double *cheby_u_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U_ZERO returns zeroes of the Chebyshev polynomial U(N)(X).
//
//  Discussion:
//
//    The I-th zero of U(N)(X) is cos((I-1)*PI/(N-1)), I = 1 to N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Output, double CHEBY_U_ZERO[N], the zeroes of U(N)(X).
//
{
  double angle;
  int i;
  static double pi = 3.141592653589793;
  double *z;

  z = new double[n];

  for ( i = 0; i < n; i++ )
  {
    angle = ( double ) ( i + 1 ) * pi / ( double ) ( n + 1 );
    z[i] = cos ( angle );
  }
  return z;
}
//****************************************************************************80

void data_to_dif ( int ntab, double xtab[], double ytab[], double diftab[] )

//****************************************************************************80
//
//  Purpose:
//
//    DATA_TO_DIF sets up a divided difference table from raw data.
//
//  Discussion:
//
//    Space can be saved by using a single array for both the DIFTAB and
//    YTAB dummy parameters.  In that case, the difference table will
//    overwrite the Y data without interfering with the computation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int NTAB, the number of pairs of points
//    (XTAB[I],YTAB[I]) which are to be used as data.
//
//    Input, double XTAB[NTAB], the X values at which data was taken.
//    These values must be distinct.
//
//    Input, double YTAB[NTAB], the corresponding Y values.
//
//    Output, double DIFTAB[NTAB], the divided difference coefficients
//    corresponding to the input (XTAB,YTAB).
//
{
  int i;
  int j;
//
//  Copy the data values into DIFTAB.
//
  for ( i = 0; i < ntab; i++ )
  {
    diftab[i] = ytab[i];
  }
//
//  Make sure the abscissas are distinct.
//
  for ( i = 0; i < ntab; i++ )
  {
    for ( j = i + 1; j < ntab; j++ )
    {
      if ( xtab[i] - xtab[j] == 0.0 )
      {
        cerr << "\n";
        cerr << "DATA_TO_DIF - Fatal error!\n";
        cerr << "  Two entries of XTAB are equal!\n";
        cerr << "  XTAB[%d] = " << xtab[i] << "\n";
        cerr << "  XTAB[%d] = " << xtab[j] << "\n";
        exit ( 1 );
      }
    }
  }
//
//  Compute the divided differences.
//
  for ( i = 1; i <= ntab - 1; i++ )
  {
    for ( j = ntab - 1; i <= j; j-- )
    {
      diftab[j] = ( diftab[j] - diftab[j-1] ) / ( xtab[j] - xtab[j-i] );
    }
  }

  return;
}
//****************************************************************************80

double *data_to_dif_new ( int ntab, double xtab[], double ytab[] )

//****************************************************************************80
//
//  Purpose:
//
//    DATA_TO_DIF_NEW sets up a divided difference table from raw data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int NTAB, the number of pairs of points
//    (XTAB[I],YTAB[I]) which are to be used as data.
//
//    Input, double XTAB[NTAB], the X values at which data was taken.
//    These values must be distinct.
//
//    Input, double YTAB[NTAB], the corresponding Y values.
//
//    Output, double DATA_TO_DIF_NEW[NTAB], the divided difference coefficients
//    corresponding to the input.
//
{
  double *diftab;
  int i;
  int j;
//
//  Make sure the abscissas are distinct.
//
  for ( i = 0; i < ntab; i++ )
  {
    for ( j = i + 1; j < ntab; j++ )
    {
      if ( xtab[i] - xtab[j] == 0.0 )
      {
        cerr << "\n";
        cerr << "DATA_TO_DIF_NEW - Fatal error!\n";
        cerr << "  Two entries of XTAB are equal!\n";
        cerr << "  XTAB[%d] = " << xtab[i] << "\n";
        cerr << "  XTAB[%d] = " << xtab[j] << "\n";
        exit ( 1 );
      }
    }
  }
//
//  Copy the Y data into DIFTAB.
//
  diftab = new double[ntab];

  for ( i = 0; i < ntab; i++ )
  {
    diftab[i] = ytab[i];
  }
//
//  Compute the divided differences.
//
  for ( i = 1; i <= ntab - 1; i++ )
  {
    for ( j = ntab - 1; i <= j; j-- )
    {
      diftab[j] = ( diftab[j] - diftab[j-1] ) 
                / ( xtab[j]   - xtab[j-i] );
    }
  }

  return diftab;
}
//****************************************************************************80

void data_to_dif_display ( int ntab, double xtab[], double ytab[],
  double diftab[] )

//****************************************************************************80
//
//  Purpose:
//
//    DATA_TO_DIF_DISPLAY sets up a divided difference table and prints intermediate data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NTAB, the number of pairs of points
//    (XTAB[I],YTAB[I]) which are to be used as data.
//
//    Input, double XTAB[NTAB], the X values at which data was taken.
//    These values must be distinct.
//
//    Input, double YTAB[NTAB], the corresponding Y values.
//
//    Output, double DIFTAB[NTAB], the divided difference coefficients
//    corresponding to the input (XTAB,YTAB).
//
{
  int i;
  int j;

  if ( !r8vec_distinct ( ntab, xtab ) )
  {
    cerr << "\n";
    cerr << "DATA_TO_DIF_DISPLAY - Fatal error!\n";
    cerr << "  Two entries of XTAB are equal!\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  The divided difference table:\n";
  cout << "\n";
  cout << "        ";
  for ( i = 0; i < ntab; i++ )
  {
    cout << setw(10) << xtab[i] << "  ";
  }
  cout << "\n";
  cout << "\n";
  cout << setw(6) << 0 << "  ";
  for ( i = 0; i < ntab; i++ )
  {
    cout << setw(10) << ytab[i] << "  ";
  }
  cout << "\n";
//
//  Copy the data values into DIFTAB.
//
  for ( i = 0; i < ntab; i++ )
  {
    diftab[i] = ytab[i];
  }
//
//  Compute the divided differences.
//
  for ( i = 1; i <= ntab - 1; i++ )
  {
    cout << setw(6) << i << "  ";
    for ( j = ntab - 1; i <= j; j-- )
    {
      diftab[j] = ( diftab[j] - diftab[j-1] ) / ( xtab[j] - xtab[j-i] );
    }
    for ( j = i; j < ntab; j++ )
    {
      cout << setw(10) << diftab[j] << "  ";
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void data_to_r8poly ( int ntab, double xtab[], double ytab[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    DATA_TO_R8POLY computes the coefficients of a polynomial interpolating data.
//
//  Discussion:
//
//    Space can be saved by using a single array for both the C and
//    YTAB parameters.  In that case, the coefficients will
//    overwrite the Y data without interfering with the computation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int NTAB, the number of data points.
//
//    Input, double XTAB[NTAB], YTAB[NTAB], the data values.
//
//    Output, double C[NTAB], the coefficients of the polynomial that passes
//    through the data (XTAB,YTAB).  C(0) is the constant term.
//
{
  if ( !r8vec_distinct ( ntab, xtab ) )
  {
    cerr << "\n";
    cerr << "DATA_TO_R8POLY - Fatal error!\n";
    cerr << "  Two entries of XTAB are equal.\n";
    exit ( 1 );
  }

  data_to_dif ( ntab, xtab, ytab, c );

  dif_to_r8poly ( ntab, xtab, c, c );

  return;
}
//****************************************************************************80

void dif_antideriv ( int ntab, double xtab[], double diftab[], int *ntab2,
  double xtab2[], double diftab2[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_ANTIDERIV integrates a polynomial in divided difference form.
//
//  Discussion:
//
//    This routine uses the divided difference representation (XTAB, DIFTAB)
//    of a polynomial to compute the divided difference representation
//    (XTAB, ANTTAB) of the antiderivative of the polynomial.
//
//    The antiderivative of a polynomial P(X) is any polynomial Q(X)
//    with the property that d/dX Q(X) = P(X).
//
//    This routine chooses the antiderivative whose constant term is zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NTAB, the size of the input table.
//
//    Input, double XTAB[NTAB], the abscissas for the divided
//    difference table.
//
//    Input, double DIFTAB[NTAB], the divided difference table.
//
//    Output, int *NTAB2, the size of the output table, which is NTAB+1.
//
//    Input, double XTAB2[NTAB2], the abscissas for the divided
//    difference table for the antiderivative.
//
//    Output, double DIFTAB2[NTAB2], the divided difference
//    table for the antiderivative.
//
{
  int i;
  double *xtab1;
  double *diftab1;
//
//  Using a temporary copy of the difference table, shift the
//  abscissas to zero.
//
  xtab1 = new double[ntab];
  diftab1 = new double[ntab];

  for ( i = 0; i < ntab; i++ )
  {
    xtab1[i] = xtab[i];
  }
  for ( i = 0; i < ntab; i++ )
  {
    diftab1[i] = diftab[i];
  }

  dif_shift_zero ( ntab, xtab1, diftab1 );

  cout << flush;
//
//  Append a final zero to XTAB.
//
  *ntab2 = ntab + 1;
  for ( i = 0; i < *ntab2; i++ )
  {
    xtab2[i] = 0.0;
  }
//
//  Get the antiderivative of the standard form polynomial.
//
  r8poly_ant_cof ( ntab, diftab1, diftab2 );

  delete [] xtab1;
  delete [] diftab1;

  return;
}
//****************************************************************************80

void dif_append ( int ntab, double xtab[], double diftab[], double xval,
  double yval, int *ntab2, double xtab2[], double diftab2[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_APPEND adds a pair of data values to a divided difference table.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NTAB, the size of the difference table.
//
//    Input, double XTAB[NTAB], the abscissas of the table.
//
//    Input, double DIFTAB[NTAB], the difference table.
//
//    Input, double XVAL, the data abscissa to be added to the table.
//
//    Input, double YVAL, the data value to be added to the table.
//
//    Output, int *NTAB2, the updated size of the difference table.
//
//    Output, double XTAB2[*NTAB2], the updated abscissas.
//
//    Output, double DIFTAB2[*NTAB2], the updated difference table.
//
{
  int i;

  *ntab2 = ntab + 1;
//
//  Move the original data up one index.
//
  for ( i = *ntab2 - 1; 1 <= i; i-- )
  {
    diftab2[i] = diftab[i-1];
    xtab2[i] = xtab[i-1];
  }
//
//  Recompute the data.
//
  xtab2[0] = xval;
  diftab2[0] = yval;

  for ( i = 1; i < *ntab2; i++ )
  {
    diftab2[i] = ( diftab2[i] - diftab2[i-1] ) / ( xtab2[i] - xtab2[0] );
  }

  return;
}
//****************************************************************************80

void dif_basis ( int ntab, double xtab[], double *diftab )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_BASIS: all Lagrange basis polynomials in divided difference form.
//
//  Discussion:
//
//    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
//    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
//    XTAB(J) for J not equal to I, and 1 when J is equal to I.
//
//    The Lagrange basis polynomials have the property that the interpolating
//    polynomial through a set of NTAB data points (XTAB,YTAB) may be
//    represented as
//
//      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
//
//    Higher order interpolation at selected points may be accomplished
//    using repeated X values, and scaled derivative values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int NTAB, the number of X data points XTAB, and the number of
//    basis polynomials to compute.
//
//    Input, double XTAB[NTAB], the X values upon which the Lagrange basis
//    polynomials are to be based.
//
//    Output, double DIFTAB[NTAB*NTAB], points to a list of NTAB * NTAB values,
//    the set of divided difference tables, stored as consecutive rows.
//    Logical row I of DIFTAB contains the table for the I-th Lagrange basis
//    polynomial.
//
{
  int i;
  int j;
  double *pointer1;
  double *pointer2;

  pointer1 = diftab;

  for ( i = 0; i < ntab; i++ )
  {
    pointer2 = pointer1;

    for ( j = 0; j < ntab; j++ )
    {
      if ( j == i )
      {
        *pointer1 = 1.0;
      }
      else
      {
        *pointer1 = 0.0;
      }
      pointer1++;
    }

    data_to_dif ( ntab, xtab, pointer2, pointer2 );
  }

  return;
}
//****************************************************************************80

void dif_basis_deriv ( int nd, double xd[], double xdp[], double ddp[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_BASIS_DERIV: Lagrange basis derivative difference tables.
//
//  Discussion:
//
//    Given ND points XD, a Lagrange basis polynomial L(J)(X) is associated
//    with each point XD(J).
//
//    This function computes a table DDP(*,*) whose J-th column contains
//    the difference table for the first derivative of L(J)(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//
//    Input, double XD[ND], the X values upon which the 
//    Lagrange basis polynomials are to be based.
//
//    Output, double XDP[ND-1], the X values upon with
//    the derivative difference table is based.  In fact, these are
//    all 0.
//
//    Output, double DDP[(ND-1)*ND], the divided difference 
//    tables for all the Lagrange basis polynomials.  Column J of DDP
//    contains the table for basis polynomial associated with XD(J).
//
{
  double *dd;
  int i;
  int j;
  double *yd;
//
//  Process the vectors one column at a time.
//
  dd = new double[nd];
  yd = new double[nd];

  for ( j = 0; j < nd; j++ )
  {
//
//  Set the data.
//
    for ( i = 0; i < nd; i++ )
    {
      yd[i] = 0.0;
    }
    yd[j] = 1.0;
//
//  Compute the divided difference table.
//
    data_to_dif ( nd, xd, yd, dd );
//
//  Compute the divided difference table for the derivative.
//
    dif_deriv_table ( nd, xd, dd, xdp, ddp + j * ( nd - 1 ) );
  }

  delete [] dd;
  delete [] yd;

  return;
}
//****************************************************************************80

void dif_basis_derivk ( int nd, double xd[], int k, double xdp[], double ddp[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_BASIS_DERIVK: Lagrange basis K-th derivative difference tables.
//
//  Discussion:
//
//    Given ND points XD, a Lagrange basis polynomial L(J)(X) is associated
//    with each point XD(J).
//
//    This function computes a table DDP(*,*) whose J-th column contains
//    the difference table for the K-th derivative of L(J)(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//
//    Input, double XD[ND], the X values upon which the 
//    Lagrange basis polynomials are to be based.
//
//    Input, int K, the index of the derivative.
//
//    Output, double XDP[ND-1], the X values upon with
//    the derivative difference table is based.  In fact, these are
//    all 0.
//
//    Output, double DDP[(ND-1)*ND], the divided difference 
//    tables for all the Lagrange basis polynomials.  Column J of DDP
//    contains the table for basis polynomial associated with XD(J).
//
{
  double *dd;
  int i;
  int j;
  double *yd;
//
//  Process the vectors one column at a time.
//
  dd = new double[nd];
  yd = new double[nd];

  for ( j = 0; j < nd; j++ )
  {
//
//  Set the data.
//
    for ( i = 0; i < nd; i++ )
    {
      yd[i] = 0.0;
    }
    yd[j] = 1.0;
//
//  Compute the divided difference table.
//
    data_to_dif ( nd, xd, yd, dd );
//
//  Compute the divided difference table for the derivative.
//
    dif_derivk_table ( nd, xd, dd, k, xdp, ddp + j * ( nd - k ) );
  }

  delete [] dd;
  delete [] yd;

  return;
}
//****************************************************************************80

void dif_basis_i ( int ival, int ntab, double xtab[], double diftab[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_BASIS_I: I-th Lagrange basis polynomial in divided difference form.
//
//  Discussion:
//
//    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
//    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
//    XTAB(J) for J not equal to I, and 1 when J is equal to I.
//
//    The Lagrange basis polynomials have the property that the interpolating
//    polynomial through a set of NTAB data points (XTAB,YTAB) may be
//    represented as
//
//      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
//
//    Higher order interpolation at selected points may be accomplished
//    using repeated X values, and scaled derivative values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int IVAL, the index of the desired Lagrange basis polynomial.
//    IVAL should be between 1 and NTAB.
//
//    Input, int NTAB, the number of data points XTAB.
//
//    Input, double XTAB[NTAB], the X values upon which the Lagrange basis
//    polynomial is to be based.
//
//    Output, double DIFTAB[NTAB], the divided difference table for the IVAL-th
//    Lagrange basis polynomial.
//
{
  int i;
//
//  Check IVAL.
//
  if ( ival < 1 || ntab < ival )
  {
    cerr << "\n";
    cerr << "DIF_BASIS_I - Fatal error!\n";
    cerr << "  IVAL must be between 1 and " << ntab << ".\n";
    cerr << "  but your value is " << ival << "\n";
    exit ( 1 );
  }
//
//  Initialize DIFTAB to Delta(I,J).
//
  for ( i = 0; i <= ntab - 1; i++ )
  {
    diftab[i] = 0.0;
  }
  diftab[ival] = 1.0;
//
//  Compute the IVAL-th Lagrange basis polynomial.
//
  data_to_dif ( ntab, xtab, diftab, diftab );

  return;
}
//****************************************************************************80

void dif_deriv_table ( int nd, double xd[], double yd[], double xdp[], 
  double ydp[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_DERIV_TABLE computes the divided difference table for a derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the size of the input table.
//
//    Input, double XD[ND], the abscissas for the divided
//    difference table.
//
//    Input, double YD[ND], the divided difference table.
//
//    Output, double XDP[ND-1], the abscissas for the divided
//    difference table for the derivative.
//
//    Output, double YDP[ND-1], the divided difference
//    table for the derivative.
//
{
  int i;
  double *xd_temp;
  double *yd_temp;
//
//  Using a temporary copy of the difference table, shift the
//  abscissas to zero.
//
  xd_temp = new double[nd];
  yd_temp = new double[nd];

  for ( i = 0; i < nd; i++ )
  {
    xd_temp[i] = xd[i];
  }
  for ( i = 0; i < nd; i++ )
  {
    yd_temp[i] = yd[i];
  }

  dif_shift_zero ( nd, xd_temp, yd_temp );
//
//  Construct the derivative.
//
  for ( i = 0; i < nd - 1; i++ )
  {
    xdp[i] = 0.0;
  }

  for ( i = 0; i < nd - 1; i++ )
  {
    ydp[i] = ( double ) ( i + 1 ) * yd_temp[i+1];
  }

  delete [] xd_temp;
  delete [] yd_temp;

  return;
}
//****************************************************************************80

void dif_derivk_table ( int nd, double xd[], double dd[], int k, 
  double xdk[], double ddk[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_DERIVK_TABLE computes the divided difference table for K-th derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the size of the input table.
//
//    Input, double XD[ND], the abscissas for the divided
//    difference table.
//
//    Input, double DD[ND], the divided difference table.
//
//    Input, int K, the index of the derivative.  0 <= K.
//
//    Input, double XDK[ND-K], the abscissas for the divided
//    difference table for the derivative.
//
//    Output, double DDK[NDP], the divided difference
//    table for the derivative.
//
{
  double *dd_temp;
  int i;
  int j;
  int ndk;
  double *xd_temp;

  if ( k < 0 )
  {
    cerr << "\n";
    cerr << "DIF_DERIVK_TABLE - Fatal error!\n";
    cerr << "  K < 0.\n";
    exit ( 1 );
  }

  if ( nd <= k )
  {
    return;
  }
//
//  Shift the abscissas to zero.
//
  ndk = nd;

  xd_temp = new double[ndk];
  dd_temp = new double[ndk];

  for ( i = 0; i < ndk; i++ )
  {
    xd_temp[i] = xd[i];
  }
  for ( i = 0; i < ndk; i++ )
  {
    dd_temp[i] = dd[i];
  }

  dif_shift_zero ( ndk, xd_temp, dd_temp );
//
//  Repeatedly differentiate.
//
  for ( j = 1; j <= k; j++ )
  {
    ndk = ndk - 1;

    for ( i = 0; i < ndk; i++ )
    {
      dd_temp[i] = ( double ) ( i + 1 ) * dd_temp[i+1];
    }
  }

  for ( i = 0; i < ndk; i++ )
  {
    ddk[i] = dd_temp[i];
    xdk[i] = 0.0;
  }

  delete [] xd_temp;
  delete [] dd_temp;

  return;
}
//****************************************************************************80

void dif_print ( int ntab, double xtab[], double diftab[], string title  )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_PRINT prints the polynomial represented by a divided difference table.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NTAB, the dimension of the arrays DIFTAB and XTAB.
//
//    Input, double XTAB[NTAB], the X values for the polynomial.
//
//    Input, double DIFTAB[NTAB], the divided difference table
//    for the polynomial.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  cout << "  p(x) =                       "
       << setw(14) << diftab[0] << "\n";

  for ( i = 1; i < ntab; i++ )
  {
    cout << "       + ( x - "
         << setw(10) << xtab[i-1] << ") * ( "
         << setw(14) << diftab[i] << "\n";
  }

  cout << "        ";
  for ( i = 1; i < ntab; i++ )
  {
    cout << ")";
  }
  cout << "\n";

  return;
}
//****************************************************************************80

void dif_shift_x ( int nd, double xd[], double yd[], double xv )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_SHIFT_X replaces one abscissa of a divided difference table.
//
//  Discussion:
//
//    This routine shifts the representation of a divided difference polynomial
//    by dropping the last X value in XD, and adding a new X value to the
//    beginning of the Xd array, suitably modifying the coefficients stored
//    in YD.
//
//    The representation of the polynomial is changed, but the polynomial itself
//    should be identical.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the number of divided difference coefficients, and
//    the number of entries in XD.
//
//    Input/output, double XD[ND], the X values used in the representation of
//    the divided difference polynomial.  After a call to this routine, the 
//    last entry of XD has been dropped, the other
//    entries have shifted up one index, and XV has been inserted at the
//    beginning of the array.
//
//    Input/output, double YD[ND], the divided difference coefficients
//    corresponding to the XD array.  On output, this array has been
//    adjusted.
//
//    Input, double XV, a new X value which is to be used in the representation
//    of the polynomial.  On output, XD[0] equals XV and the representation
//    of the polynomial has been suitably changed.
//    Note that XV does not have to be distinct from any of the original XD
//    values.
//
{
  int i;
//
//  Recompute the divided difference coefficients.
//
  for ( i = nd - 2; 0 <= i; i-- )
  {
    yd[i] = yd[i] + ( xv - xd[i] ) * yd[i+1];
  }
//
//  Shift the X values up one position and insert XV.
//
  for ( i = nd - 1; 0 < i; i-- )
  {
    xd[i] = xd[i-1];
  }

  xd[0] = xv;

  return;
}
//****************************************************************************80

void dif_shift_zero ( int nd, double xd[], double yd[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_SHIFT_ZERO shifts a divided difference table so all abscissas are zero.
//
//  Discussion:
//
//    When the abscissas are changed, the coefficients naturally
//    must also be changed.
//
//    The resulting pair (XD, YD) still represents the
//    same polynomial, but the entries in YD are now the
//    standard polynomial coefficients.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the length of the XD and YD arrays.
//
//    Input/output, double XD[ND], the X values that correspond to the
//    divided difference table.  On output, XD contains only zeroes.
//
//    Input/output, double YD[ND], the divided difference table
//    for the polynomial.  On output, YD is also
//    the coefficient array for the standard representation
//    of the polynomial.
//
{
  int i;
  int j;

  for ( j = 1; j <= nd; j++ )
  {
//
//  Recompute the divided difference coefficients.
//
    for ( i = nd - 2; 0 <= i; i-- )
    {
      yd[i] = yd[i] - xd[i] * yd[i+1];
    }
//
//  Shift the XD values up one position and insert XV.
//
    for ( i = nd - 1; 0 < i; i-- )
    {
      xd[i] = xd[i-1];
    }
    xd[0] = 0.0;
  }

  return;
}
//****************************************************************************80

void dif_to_r8poly ( int ntab, double xtab[], double diftab[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_TO_R8POLY converts a divided difference table to a standard polynomial.
//
//  Discussion:
//
//    The vector DIFTAB, containing the divided difference polynomial
//    coefficients is overwritten with the standard form polynomial
//    coefficients, but the abscissa vector XTAB is unchanged.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int NTAB, the number of coefficients, and abscissas.
//
//    Input, double XTAB[NTAB], the X values used in the divided difference
//    representation of the polynomial.
//
//    Input, double DIFTAB[NTAB], the divided difference table.
//
//    Output, double C[NTAB], the standard form polyomial coefficients.
//    C[0] is the constant term, and C[NTAB-1] is the coefficient
//    of X**(NTAB-1).
//
{
  int i;
  int j;

  for ( i = 0; i < ntab; i++ )
  {
    c[i] = diftab[i];
  }
//
//  Recompute the divided difference coefficients.
//
  for ( j = 1; j <= ntab - 1; j++ )
  {
    for ( i = 1; i <= ntab - j; i++ )
    {
      c[ntab-i-1] = c[ntab-i-1] - xtab[ntab-i-j] * c[ntab-i];
    }
  }

  return;
}
//****************************************************************************80

double dif_val ( int ntab, double xtab[], double diftab[], double xv )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_VAL evaluates a divided difference polynomial at a point.
//
//  Discussion:
//
//    DATA_TO_DIF must be called first to set up the divided difference table.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, integer NTAB, the number of divided difference
//    coefficients in DIFTAB, and the number of points XTAB.
//
//    Input, double XTAB[NTAB], the X values upon which the
//    divided difference polynomial is based.
//
//    Input, double DIFTAB[NTAB], the divided difference table.
//
//    Input, double XV, a value of X at which the polynomial
//    is to be evaluated.
//
//    Output, double DIF_VAL, the value of the polynomial at XV.
//
{
  int i;
  double yv;

  yv = diftab[ntab-1];
  for ( i = 2; i <= ntab; i++ )
  {
    yv = diftab[ntab-i] + ( xv - xtab[ntab-i] ) * yv;
  }

  return yv;
}
//****************************************************************************80

double *dif_vals ( int nd, double xd[], double yd[], int nv, double xv[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIF_VALS evaluates a divided difference polynomial at a set of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer, 2001,
//    ISBN: 0387953663,
//    LC: QA1.A647.v27.
//
//  Parameters:
//
//    Input, int ND, the order of the difference table.
//
//    Input, double XD[ND], the X values of the difference table.
//
//    Input, double YD[ND], the divided differences.
//
//    Input, int NV, the number of evaluation points.
//
//    Input, double XV[NV], the evaluation points.
//
//    Output, double DIF_VALS[NV], the value of the divided difference
//    polynomial at the evaluation points.
//
{
  int i;
  int j;
  double *yv;

  yv = new double[nv];

  for ( j = 0; j < nv; j++ )
  {
    yv[j] = yd[nd-1];
    for ( i = 2; i <= nd; i++ )
    {
      yv[j] = yd[nd-i] + ( xv[j] - xd[nd-i] ) * yv[j];
    }
  }
  return yv;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

double *lagrange_rule ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_RULE computes the weights of a Lagrange interpolation rule.
//
//  Discussion:
//
//    Given N abscissas X, an arbitrary function F(X) can be
//    interpolated by a polynomial P(X) of order N (and degree N-1)
//    using weights that depend only on X.
//
//    Standard Lagrange interpolation can be rewritten into this form,
//    which is more economical than evaluating the individual Lagrange
//    basis polynomials.
//
//    If we define
//
//      W(I) = 1 / product ( 1 <= J <= N, J /= I ) ( X(J) - X(I) )
//
//    then
//
//      P(XV) = sum ( 1 <= I <= N ) W(I) * F( X(I) ) / ( XV - X(I) )
//            / sum ( 1 <= I <= N ) W(I)             / ( XV - X(I) )
//
//    except when XV = X(J), for some J, when we set:
//
//      P(X(J)) = F(X(J))
//
//  Modified:
//
//    24 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jean-Paul Berrut, Lloyd Trefethen,
//    Barycentric Lagrange Interpolation,
//    SIAM Review,
//    Volume 46, Number 3, September 2004, pages 501-517.
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Input, double X[N], the abscissas of the rule.
//
//    Output, double LAGRANGE_RULE[N], the weights of the rule.
//
{
  int i;
  int j;
  double *w;

  w = new double[n];

  for ( i = 0; i < n; i++ )
  {
    w[i] = 1.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i != j )
      {
        w[j] = w[j] / ( x[i] - x[j] );
      }
    }
  }
  return w;
}
//****************************************************************************80

double lagrange_sum ( int n, double x[], double w[], double y[], double xv )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_SUM carries out a Lagrange interpolation rule.
//
//  Discussion:
//
//    It is assumed that LAGRANGE_RULE has already been called to compute
//    the appropriate weights for the given set of abscissas.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2001
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jean-Paul Berrut, Lloyd Trefethen,
//    Barycentric Lagrange Interpolation,
//    SIAM Review,
//    Volume 46, Number 3, September 2004, pages 501-517.
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Input, double X[N], the abscissas of the rule.
//
//    Input, double W[N], the weights of the rule.
//
//    Input, double Y[N], the function values at the abscissas.
//
//    Input, double XV, a point where an interpolated value is
//    needed.
//
//    Output, double LAGRANGE_SUM, the interpolated function value.
//
{
  double bot;
  int  i;
  double top;
  double yv;

  for ( i = 0; i < n; i++ )
  {
    if ( xv == x[i] )
    {
      yv = y[i];
      return yv;
    }
  }

  top = 0.0;
  bot = 0.0;
  
  for ( i = 0; i < n; i++ )
  {
    top = top + w[i] * y[i] / ( xv - x[i] );
    bot = bot + w[i]        / ( xv - x[i] );
  }

  yv = top / bot;

  return yv;
}
//****************************************************************************80

double lagrange_val ( int n, double x[], double y[], double xv )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_VAL applies a naive form of Lagrange interpolation.
//
//  Discussion:
//
//    Given N abscissas X, an arbitrary function Y(X) can be
//    interpolated by a polynomial P(X) of order N (and degree N-1)
//    using Lagrange basis polynomials of degree N-1.
//
//    Standard Lagrange interpolation can be rewritten into this form,
//    which is more economical than evaluating the individual Lagrange
//    basis polynomials.
//
//    If we define
//
//      L(I)(XV) = product ( 1 <= J <= N, J /= I )
//        ( XV - X(J) ) / ( X(I) - X(J) )
//
//    then
//
//      P(XV) = sum ( 1 <= I <= N ) Y( X(I) ) * L(I)(XV)
//
//    Applying this form of the interpolation rule directly involves 
//    about N^2 work.  There are more efficient forms of the rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data points.
//
//    Input, double X[N], the abscissas.
//
//    Input, double Y[N], the function values at the abscissas.
//
//    Input, double XV, a point where an interpolated value is
//    needed.
//
//    Output, double LAGRANGE_VAL, the interpolated function value.
//
{
  int i;
  int j;
  double poly;
  double yv;

  yv = 0.0;

  for ( i = 0; i < n; i++ )
  {
    poly = 1.0;
    for ( j = 0; j < n; j++ )
    {
      if ( j != i )
      {
        poly = poly * ( xv - x[j] ) / ( x[i] - x[j] );
      }
    }
    yv = yv + y[i] * poly;
  }
  return yv;
}
//****************************************************************************80

void nc_rule ( int norder, double a, double b, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    NC_RULE computes the weights of a Newton-Cotes quadrature rule.
//
//  Discussion:
//
//    For the interval [A,B], the Newton-Cotes quadrature rule estimates
//
//      Integral ( A <= X <= B ) F(X) dX
//
//    using NORDER equally spaced abscissas XTAB(I) and a weight vector
//    WEIGHT(I):
//
//      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
//
//    For the CLOSED rule, the abscissas include the points A and B.
//    For the OPEN rule, the abscissas do not include A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NORDER, the order of the rule.
//
//    Input, double A, B, the left and right endpoints of the interval
//    over which the quadrature rule is to be applied.
//
//    Input, double XTAB[NORDER], the abscissas of the rule.
//
//    Output, double WEIGHT[NORDER], the weights of the rule.
//
{
  int i;
  double *poly_cof;
  double yvala;
  double yvalb;
//
//  Allocate temporary space for POLY_COF.
//
  poly_cof = new double[norder];

  for ( i = 1; i <= norder; i++ )
  {
//
//  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
//  and zero at the other nodes.
//
    r8poly_basis_1 ( i, norder, xtab, poly_cof );
//
//  Evaluate the antiderivative of the polynomial at the left and
//  right endpoints.
//
    yvala = r8poly_ant_val ( norder-1, poly_cof, a );

    yvalb = r8poly_ant_val ( norder-1, poly_cof, b );

    weight[i-1] = yvalb - yvala;
  }
//
//  Free up POLY_COF space.
//
  delete [] poly_cof;

  return;
}
//****************************************************************************80

void ncc_rule ( int norder, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCC_RULE computes the coefficients of a Newton-Cotes closed quadrature rule.
//
//  Discussion:
//
//    For the interval [-1,1], the Newton-Cotes quadrature rule estimates
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    using NORDER equally spaced abscissas XTAB(I) and a weight vector
//    WEIGHT(I):
//
//      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
//
//    For the CLOSED rule, the abscissas include A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NORDER, the order of the rule.
//
//    Output, double XTAB[NORDER], the abscissas of the rule.
//
//    Output, double WEIGHT[NORDER], the weights of the rule.
//
{
  double a;
  double b;
  int i;
//
//  Compute a closed quadrature rule.
//
  a = -1.0;
  b = 1.0;

  for ( i = 1; i <= norder; i++ )
  {
    xtab[i-1] = ( ( double ) ( norder - i ) * a + ( double ) ( i - 1 ) * b )
      / ( double ) ( norder - 1 );
  }

  nc_rule ( norder, a, b, xtab, weight );

  return;
}
//****************************************************************************80

void nco_rule ( int norder, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_RULE computes the coefficients of a Newton-Cotes open quadrature rule.
//
//  Discussion:
//
//    For the interval [-1,1], the Newton-Cotes quadrature rule estimates
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    using NORDER equally spaced abscissas XTAB(I) and a weight vector
//    WEIGHT(I):
//
//      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
//
//    For the OPEN rule, the abscissas do not include A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NORDER, the order of the rule.
//
//    Output, double XTAB[NORDER], the abscissas of the rule.
//
//    Output, double WEIGHT[NORDER], the weights of the  rule.
//
{
  double a;
  double b;
  int i;

  a = -1.0;
  b = 1.0;

  for ( i = 1; i <= norder; i++ )
  {
    xtab[i-1] = ( ( double ) ( norder + 1 - i ) * a + ( double ) ( i ) * b )
      / ( double ) ( norder + 1 );
  }

  nc_rule ( norder, a, b, xtab, weight );

  return;
}
//****************************************************************************80

void r8_swap ( double *x, double *y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SWAP swaps two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, double *X, *Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  double z;

  z = *x;
  *x = *y;
  *y = z;

  return;
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i - 1 << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j - 1 << ":";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8poly_ant_cof ( int n, double poly_cof[], double poly_cof2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_ANT_COF integrates an R8POLY in standard form.
//
//  Discussion:
//
//    The antiderivative of a polynomial P(X) is any polynomial Q(X)
//    with the property that d/dX Q(X) = P(X).
//
//    This routine chooses the antiderivative whose constant term is zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double POLY_COF[N], the polynomial coefficients.
//    POLY_COF[0] is the constant term, and POLY_COF[N-1] is the
//    coefficient of X**(N-1).
//
//    Output, double POLY_COF2[N+1], the coefficients of the antiderivative
//    polynomial, in standard form.  The constant term is set to zero.
//
{
  int i;
//
//  Set the constant term.
//
  poly_cof2[0] = 0.0;
//
//  Integrate the polynomial.
//
  for ( i = 1; i <= n; i++ )
  {
    poly_cof2[i] = poly_cof[i-1] / ( double ) i;
  }

  return;
}
//****************************************************************************80

double r8poly_ant_val ( int n, double poly_cof[], double xval )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_ANT_VAL evaluates the antiderivative of an R8POLY in standard form.
//
//  Discussion:
//
//    The constant term of the antiderivative is taken to be zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double POLY_COF[N], the polynomial coefficients.  POLY_COF[0]
//    is the constant term, and POLY_COF[N-1] is the coefficient of X**(N-1).
//
//    Input, double XVAL, the point where the antiderivative is to be
//    evaluated.
//
//    Output, double R8POLY_ANT_VAL, the value of the antiderivative of the polynomial
//    at XVAL.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    value = ( value + poly_cof[i] / ( double ) ( i + 1 ) ) * xval;
  }

  return value;
}
//****************************************************************************80

void r8poly_basis ( int ntab, double xtab[], double *poly_cof )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_BASIS computes all Lagrange basis polynomials in standard form.
//
//  Discussion:
//
//    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
//    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
//    XTAB(J) for J not equal to I, and 1 when J is equal to I.
//
//    The Lagrange basis polynomials have the property that the interpolating
//    polynomial through a set of NTAB data points (XTAB,YTAB) may be
//    represented as
//
//      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
//
//    Higher order interpolation at selected points may be accomplished
//    using repeated X values, and scaled derivative values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NTAB, the number of data points XTAB.
//
//    Input, double XTAB[NTAB], the X values upon which the Lagrange basis
//    polynomial is to be based.
//
//    Output, double *POLY_COF, points to NTAB * NTAB values.  The polynomial
//    coefficients for the I-th Lagrange basis polynomial are stored in
//    (logical) row I.  POLY_COF[0,*] is the constant term, and POLY_COF[NTAB-1,*] is
//    the coefficient of X**(NTAB-1).
//
{
  int i;
  int j;
  double *pointer1;
  double *pointer2;

  pointer1 = poly_cof;

  for ( i = 0; i < ntab; i++ )
  {
    pointer2 = pointer1;

    for ( j = 0; j < ntab; j++ )
    {
      if ( j == i )
      {
        *pointer1 = 1.0;
      }
      else
      {
        *pointer1 = 0.0;
      }
      pointer1++;
    }
//
//  Compute the divided difference table for the IVAL-th Lagrange basis
//  polynomial.
//
    data_to_dif ( ntab, xtab, pointer2, pointer2 );
//
//  Convert the divided difference table coefficients to standard polynomial
//  coefficients.
//
    dif_to_r8poly ( ntab, xtab, pointer2, pointer2 );
  }

  return;
}
//****************************************************************************80

void r8poly_basis_1 ( int ival, int ntab, double xtab[], double poly_cof[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_BASIS_1 computes the I-th Lagrange basis polynomial in standard form.
//
//  Discussion:
//
//    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
//    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
//    XTAB(J) for J not equal to I, and 1 when J is equal to I.
//
//    The Lagrange basis polynomials have the property that the interpolating
//    polynomial through a set of NTAB data points (XTAB,YTAB) may be
//    represented as
//
//      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
//
//    Higher order interpolation at selected points may be accomplished
//    using repeated X values, and scaled derivative values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, the index of the desired Lagrange basis polynomial.
//    IVAL should be between 1 and NTAB.
//
//    Input, int NTAB, the number of data points XTAB.
//
//    Input, double XTAB[NTAB], the X values upon which the Lagrange basis
//    polynomial is to be based.
//
//    Output, double POLY_COF[NTAB], the polynomial coefficients for the
//    IVAL-th Lagrange basis polynomial.
//
{
  int i;
//
//  Check IVAL.
//
  if ( ival < 1 || ntab < ival )
  {
    cerr << "\n";
    cerr << "R8POLY_BASIS_1 - Fatal error!\n";
    cerr << "  IVAL must be between 1 and " << ntab << ".\n";
    cerr << "  but your value is " << ival << ".\n";
    exit ( 1 );
  }
//
//  Initialize POLY_COF to the IVAL-th column of the identity matrix.
//
  for ( i = 0; i <= ntab - 1; i++ )
  {
    poly_cof[i] = 0.0;
  }
  poly_cof[ival-1] = 1.0;
//
//  Compute the divided difference table for the IVAL-th Lagrange basis
//  polynomial.
//
  data_to_dif ( ntab, xtab, poly_cof, poly_cof );
//
//  Convert the divided difference table coefficients to standard polynomial
//  coefficients.
//
  dif_to_r8poly ( ntab, xtab, poly_cof, poly_cof );

  return;
}
//****************************************************************************80

void r8poly_der_cof ( int n, double poly_cof[], double poly_cof2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_DER_COF computes the coefficients of the derivative of a real polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double POLY_COF[N], the coefficients of the polynomial to
//    be differentiated.  POLY_COF[0] is the constant term, and
//    POLY_COF[N-1] is the coefficient of X**(N-1).
//
//    Output, double POLY_COF2[N-1], the coefficients of the derivative of
//    the polynomial.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    poly_cof2[i] = ( double ) ( i + 1 ) * poly_cof[i+1];
  }

  return;
}
//****************************************************************************80

double r8poly_der_val ( int n, double poly_cof[], double xval )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_DER_VAL evaluates the derivative of a real polynomial in standard form.
//
//  Discussion:
//
//    A polynomial in standard form, with coefficients POLY_COF(*),
//    may be written:
//
//    P(X) = POLY_COF[0]
//         + POLY_COF[1] * X
//         ...
//         + POLY_COF[N-1] * X**(N-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double POLY_COF[N], the polynomial coefficients.  POLY_COF[0]
//    is the constant term, and POLY_COF[N-1] is the coefficient of
//    X**(N-1).
//
//    Input, double XVAL, a value where the derivative of the polynomial
//    is to be evaluated.
//
//    Output, double R8POLY_DER_VAL, the value of the derivative of the polynomial
//    at XVAL.
//
{
  int i;
  double value;

  value = ( double ) ( n - 1 ) * poly_cof[n-1];

  for ( i = n - 2; 1 <= i; i-- )
  {
    value = value * xval + ( double ) i * poly_cof[i];
  }

  return value;
}
//****************************************************************************80

int r8poly_order ( int na, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_ORDER returns the order of a polynomial.
//
//  Discussion:
//
//    The order of a polynomial is the degree plus 1.
//
//    The order of a constant polynomial is 1.  The order of the
//    zero polynomial is debatable, but this routine returns the
//    order as 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the order of the polynomial (ignoring zero coefficients).
//
//    Input, double A[NA], the coefficients of the polynomial.
//
//    Output, int R8POLY_ORDER, the degree of the polynomial.
//
{
  int value;

  value = na;

  while ( 1 < value )
  {
    if ( a[value-1] != 0.0 )
    {
      return value;
    }
    value = value - 1;
  }

  return value;
}
//****************************************************************************80

void r8poly_print ( int n, double poly_cof[], string title  )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_PRINT prints out a real polynomial in standard form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double POLY_COF[N], the polynomial coefficients.
//    POLY_COF[0] is the constant term and
//    POLY_COF[N-1] is the coefficient of X^(N-1).
//
//    Input, string TITLE, a title.
//
{
  int i;
  int n2;

  cout << "\n";
  cout << title << "\n";
//
//  Determine the "real" degree of the polynomial.
//  N2 is the index of the highest order nonzero coefficient.
//
  n2 = r8poly_order ( n, poly_cof );

  cout << "\n";
  cout << "  p(x) = ";

  if ( poly_cof[n2-1] < 0.0 )
  {
    cout << "- ";
  }
  else
  {
    cout << "  ";
  }

  cout << setw(8) << fabs ( poly_cof[n2-1] );

  if ( 2 <= n2 )
  {
    cout << " * x";
  }

  if ( 3 <= n2 )
  {
    cout << " ^ " << n2 - 1;
  }
  cout << "\n";

  for ( i = n2 - 1; 1 <= i; i-- )
  {

    if ( poly_cof[i-1] < 0.0 )
    {
      cout << "       - ";
    }
    else if ( poly_cof[i-1] == 0.0 )
    {
      break;
    }
    else if ( 0.0 <= poly_cof[i-1] )
    {
      cout << "       + ";
    }

    cout << setw(10) << fabs ( poly_cof[i-1] );

    if ( 2 <= i )
    {
      cout <<" * x";
    }

    if ( 3 <= i )
    {
      cout << " ^ " << i - 1;
    }

    cout << "\n";

  }

  return;
}
//****************************************************************************80

void r8poly_shift ( double scale, double shift, int n, double poly_cof[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_SHIFT adjusts the coefficients of a polynomial for a new argument.
//
//  Discussion:
//
//    Assuming P(X) is a polynomial in the argument X, of the form:
//
//      P(X) =
//          C(N-1) * X^(N-1)
//        + ...
//        + C(1) * X
//        + C(0),
//
//    and that Z is related to X by the formula:
//
//      Z = SCALE * X + SHIFT
//
//    then this routine computes coefficients C for the polynomial Q(Z):
//
//      Q(Z) =
//          C(N-1) * Z^(N-1)
//        + ...
//        + C(1) * Z
//        + C(0)
//
//    so that:
//
//      Q(Z(X)) = P(X)
//
//  Example:
//
//    P(X) = 2 * X^2 - X + 6
//
//    Z = 2.0 * X + 3.0
//
//    Q(Z) = 0.5 *         Z^2 -  3.5 * Z + 12
//
//    Q(Z(X)) = 0.5 * ( 4.0 * X^2 + 12.0 * X +  9 )
//            - 3.5 * (               2.0 * X +  3 )
//                                            + 12
//
//            = 2.0         * X^2 -  1.0 * X +  6
//
//            = P(X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double SHIFT, SCALE, the shift and scale applied to X,
//    so that Z = SCALE * X + SHIFT.
//
//    Input, int N, the order of the polynomial.
//
//    Input/output, double POLY_COF[N].
//    On input, the coefficient array in terms of the X variable.
//    On output, the coefficient array in terms of the Z variable.
//
{
  int i;
  int j;

  for ( i = 1; i <= n; i++ )
  {
    for ( j = i+1; j <= n; j++ )
    {
      poly_cof[j-1] = poly_cof[j-1] / scale;
    }
  }

  for ( i = 1; i <= n; i++ )
  {
    for ( j = n-1; i <= j; j-- )
    {
      poly_cof[j-1] = poly_cof[j-1] - shift * poly_cof[j];
    }
  }

  return;
}
//****************************************************************************80

double r8poly_val_horner ( int n, double poly_cof[], double xval )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_VAL_HORNER evaluates a real polynomial in standard form.
//
//  Discussion:
//
//    A polynomial in standard form, with coefficients POLY_COF(*),
//    may be written:
//
//    P(X) = POLY_COF[0]
//         + POLY_COF[1] * X
//         ...
//         + POLY_COF[N-1] * X^(N-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double POLY_COF[N], the polynomial coefficients.  POLY_COF[0]
//    is the constant term, and POLY_COF[N-1] is the coefficient of
//    X^(N-1).
//
//    Input, double XVAL, a value where the polynomial is to be evaluated.
//
//    Output, double R8POLY_VAL_HORNER, the value of the polynomial at XVAL.
//
{
  int i;
  double value;

  value = poly_cof[n-1];

  for ( i = n - 2; 0 <= i; i-- )
  {
    value = value * xval + poly_cof[i];
  }

  return value;
}
//****************************************************************************80

bool r8vec_distinct ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[N], the vector to be checked.
//
//    Output, bool R8VEC_DISTINCT is true if all N elements of X are distinct.
//
{
  int i;
  int j;

  for ( i = 1; i <= n - 1; i++ )
  {
    for ( j = 1; j <= i - 1; j++ )
    {
      if ( x[i] == x[j] )
      {
        return false;
      }
    }
  }

  return true;
}
//****************************************************************************80

void r8vec_indicator ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDICATOR sets an R8VEC to the indicator vector {1,2,3...}.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, double A[N], the array to be initialized.
//
{
  int i;

  for ( i = 0; i <= n - 1; i++ )
  {
    a[i] = ( double ) ( i + 1 );
  }

  return;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 December 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n - 1; i++ )
  {
    cout << setw(6)  << i + 1  << "  "
         << setw(10) << a[i] << "\n";
  }

  return;
}
//****************************************************************************80

void roots_to_dif ( int nroots, double roots[], int *ntab, double xtab[],
  double diftab[] )

//****************************************************************************80
//
//  Purpose:
//
//    ROOTS_TO_DIF sets a divided difference table for a polynomial from its roots.
//
//  Discussion:
//
//    This turns out to be a simple task, because of two facts:
//
//    * The divided difference polynomial of one smaller degree which
//      passes through the values ( ROOT(I), 0 ) is the zero polynomial,
//      and hence has a zero divided difference table.
//
//    * We want a polynomial of one degree higher, but we don't want it
//      to pass through an addditional point.  Instead, we specify that
//      the polynomial is MONIC.  This means that the divided difference
//      table is almost the same as for the zero polynomial, except that
//      there is one more pair of entries, an arbitrary X value, and
//      a Y value of 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NROOTS, is the number of roots.
//
//    Input, double ROOTS[NROOTS], the roots of the polynomial.
//
//    Output, int *NTAB, the size of the divided difference table.
//    This is NROOTS+1
//
//    Output, double XTAB[NTAB], the abscissas of the table.
//
//    Output, double DIFTAB[NTAB], the divided difference table.
//
{
  int i;

  *ntab = nroots + 1;
//
//  Build the appropriate difference table for the polynomial
//  through ( ROOTS(I), 0 ) of degree NTAB-2.
//
  for ( i = 0; i < *ntab - 1; i++ )
  {
    diftab[i] = 0.0;
  }
  for ( i = 0; i < *ntab - 1; i++ )
  {
    xtab[i] = roots[i];
  }
//
//  Append the extra data to make a monic polynomial of degree NTAB-1
//  which is zero at the NTAB-1 roots.
//
  xtab[*ntab-1] = 0.0;
  diftab[*ntab-1] = 1.0;

  return;
}
//****************************************************************************80

void roots_to_r8poly ( int nroots, double roots[], int *nc, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NROOTS, the number of roots specified.
//
//    Input, double ROOTS[NROOTS], the roots.
//
//    Output, integer *NC, the order of the polynomial, which will be NROOTS+1.
//
//    Output, double C[*NC], the coefficients of the polynomial.
//
{
  int i;
  int j;
  double *xtab;

  *nc = nroots + 1;
//
//  Initialize C to (0, 0, ..., 0, 1).
//  Essentially, we are setting up a divided difference table.
//
  xtab = new double[nroots+1];
  for ( i = 0; i < *nc-1; i++ )
  {
    xtab[i] = roots[i];
  }
  xtab[*nc-1] = 0.0;

  for ( i = 0; i < *nc - 1; i++ )
  {
    c[i] = 0.0;
  }
  c[*nc-1] = 1.0;
//
//  Convert to standard polynomial form by shifting the abscissas
//  of the divided difference table to 0.
//
  dif_shift_zero ( *nc, xtab, c );

  delete [] xtab;

  return;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
