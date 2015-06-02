# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "toms886.hpp"

//****************************************************************************80

void cheb ( int deg, double pt, double tcheb[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEB computes normalized Chebyshev polynomials.
//
//  Discussion:
//
//    This subroutine computes the array TCHEB of normalized Chebyshev 
//    polynomials from degree 0 to DEG:
//      T_0(x)=1, 
//      T_j(x) = sqrt(2) * cos ( j * acos(x) ) 
//    at the point x = PT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Parameters:
//
//    Input, int DEG, the degree.
//    0 <= DEG.
//
//    Input, double PT, the evaluation point.
//
//    Output, double TCHEB[DEG+1], the value of the normalized
//    Chebyshev polynomials of degrees 0 through DEG at the point PT.
//
{
  int j;
  const double sqrt2 = 1.4142135623730951;

  if ( deg < 0 )
  {
    return;
  }

  tcheb[0] = 1.0;

  if ( deg < 1 )
  {
    return;
  }

  tcheb[1] = sqrt2 * pt;
 
  if ( deg < 2 )
  {
    return;
  }

  tcheb[2] = 2.0 * pt * tcheb[1] - sqrt2 * tcheb[0];
//
//  Chebyshev recurrence.
//
  for ( j = 3; j <= deg; j++ )
  {
    tcheb[j] = 2.0 * pt * tcheb[j-1] - tcheb[j-2];
  }

  return;
}
//****************************************************************************80

void dgemm ( char transa, char transb, int m, int n, int k, 
  double alpha, double a[], int lda, double b[], int ldb, double beta, 
  double c[], int ldc )

//****************************************************************************80
//
//  Purpose:
//
//    DGEMM computes C = alpha * A * B and related operations.
//
//  Discussion:
//
//    DGEMM performs one of the matrix-matrix operations
//
//     C := alpha * op ( A ) * op ( B ) + beta * C,
//
//    where op ( X ) is one of
//
//      op ( X ) = X   or   op ( X ) = X',
//
//    ALPHA and BETA are scalars, and A, B and C are matrices, with op ( A )
//    an M by K matrix, op ( B ) a K by N matrix and C an N by N matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, char TRANSA, specifies the form of op( A ) to be used in
//    the matrix multiplication as follows:
//    'N' or 'n', op ( A ) = A.
//    'T' or 't', op ( A ) = A'.
//    'C' or 'c', op ( A ) = A'.
//
//    Input, char TRANSB, specifies the form of op ( B ) to be used in
//    the matrix multiplication as follows:
//    'N' or 'n', op ( B ) = B.
//    'T' or 't', op ( B ) = B'.
//    'C' or 'c', op ( B ) = B'.
//
//    Input, int M, the number of rows of the  matrix op ( A ) and of the  
//    matrix C.  0 <= M.
//
//    Input, int N, the number  of columns of the matrix op ( B ) and the 
//    number of columns of the matrix C.  0 <= N.
//
//    Input, int K, the number of columns of the matrix op ( A ) and the 
//    number of rows of the matrix op ( B ).  0 <= K.
//
//    Input, double ALPHA, the scalar multiplier 
//    for op ( A ) * op ( B ).
//
//    Input, double A(LDA,KA), where:
//    if TRANSA is 'N' or 'n', KA is equal to K, and the leading M by K
//    part of the array contains A;
//    if TRANSA is not 'N' or 'n', then KA is equal to M, and the leading
//    K by M part of the array must contain the matrix A.
//
//    Input, int LDA, the first dimension of A as declared in the calling 
//    routine.  When TRANSA = 'N' or 'n' then LDA must be at least max ( 1, M ), 
//    otherwise LDA must be at least max ( 1, K ).
//
//    Input, double B(LDB,KB), where:
//    if TRANSB is 'N' or 'n', kB is N, and the leading K by N 
//    part of the array contains B;
//    if TRANSB is not 'N' or 'n', then KB is equal to K, and the leading
//    N by K part of the array must contain the matrix B.
//
//    Input, int LDB, the first dimension of B as declared in the calling 
//    routine.  When TRANSB = 'N' or 'n' then LDB must be at least max ( 1, K ), 
//    otherwise LDB must be at least max ( 1, N ).
//
//    Input, double BETA, the scalar multiplier for C.
//
//    Input, double C(LDC,N).
//    Before entry, the leading M by N part of this array must contain the 
//    matrix C, except when BETA is 0.0, in which case C need not be set 
//    on entry.
//    On exit, the array C is overwritten by the M by N matrix
//      alpha * op ( A ) * op ( B ) + beta * C.
//
//    Input, int LDC, the first dimension of C as declared in the calling 
//    routine.  max ( 1, M ) <= LDC.
//
{
  int i;
  int info;
  int j;
  int l;
  int ncola;
  int nrowa;
  int nrowb;
  bool nota;
  bool notb;
  double temp;
//
//  Set NOTA and NOTB as true if A and B respectively are not
//  transposed and set NROWA, NCOLA and NROWB as the number of rows
//  and columns of A and the number of rows of B respectively.
//
  nota = ( ( transa == 'N' ) || ( transa == 'n' ) );

  if ( nota )
  {
    nrowa = m;
    ncola = k;
  }
  else
  {
    nrowa = k;
    ncola = m;
  }

  notb = ( ( transb == 'N' ) || ( transb == 'n' ) );

  if ( notb )
  {
    nrowb = k;
  }
  else
  {
    nrowb = n;
  }
//
//  Test the input parameters.
//
  info = 0;

  if ( ! ( transa == 'N' || transa == 'n' ||
           transa == 'C' || transa == 'c' ||
           transa == 'T' || transa == 't' ) )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input TRANSA has an illegal value.\n";
    exit ( 1 );
  }

  if ( ! ( transb == 'N' || transb == 'n' ||
           transb == 'C' || transb == 'c' ||
           transb == 'T' || transb == 't' ) )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input TRANSB has an illegal value.\n";
    exit ( 1 );
  }

  if ( m < 0 )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input M has an illegal value.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input N has an illegal value.\n";
    exit ( 1 );
  }

  if ( k  < 0 )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input K has an illegal value.\n";
    exit ( 1 );
  }

  if ( lda < i4_max ( 1, nrowa ) )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input LDA has an illegal value.\n";
    exit ( 1 );
  }

  if ( ldb < i4_max ( 1, nrowb ) )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input LDB has an illegal value.\n";
    exit ( 1 );
  }

  if ( ldc < i4_max ( 1, m ) )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input LDC has an illegal value.\n";
    exit ( 1 );
  }
//
//  Quick return if possible.
//
  if ( m == 0 )
  {
    return;
  }

  if ( n == 0 )
  {
    return;
  }

  if ( ( alpha == 0.0 || k == 0 ) && ( beta == 1.0 ) )
  {
    return;
  }
//
//  And if alpha is 0.0.
//
  if ( alpha == 0.0 )
  {
    if ( beta == 0.0 )
    {
      for ( j = 0; j < n; j++ )
      {
        for ( i = 0; i < m; i++ )
        {
          c[i+j*ldc] = 0.0;
        }
      }
    }
    else
    {
      for ( j = 0; j < n; j++ )
      {
        for ( i = 0; i < m; i++ )
        {
          c[i+j*ldc] = beta * c[i+j*ldc];
        }
      }
    }
    return;
  }
//
//  Start the operations.
//
  if ( notb )
  {
//
//  Form  C := alpha*A*B + beta*C.
//
    if ( nota )
    {
      for ( j = 0; j < n; j++ )
      {
        if ( beta == 0.0 )
        {
          for ( i = 0; i < m; i++ )
          {
            c[i+j*ldc] = 0.0;
          }
        }
        else if ( beta != 1.0 )
        {
          for ( i = 0; i < m; i++ )
          {
            c[i+j*ldc] = beta * c[i+j*ldc];
          }
        }

        for ( l = 0; l < k; l++ )
        {
          if ( b[l+j*ldb] != 0.0 )
          {
            temp = alpha * b[l+j*ldb];
            for ( i = 0; i < m; i++ )
            {
              c[i+j*ldc] = c[i+j*ldc] + temp * a[i+l*lda];
            }
          }
        }

      }
    }
//
//  Form  C := alpha*A'*B + beta*C
//
    else
    {
      for ( j = 0; j < n; j++ )
      {
        for ( i = 0; i < m; i++ )
        {
          temp = 0.0;
          for ( l = 0; l < k; l++ )
          {
            temp = temp + a[l+i*lda] * b[l+j*ldb];
          }

          if ( beta == 0.0 )
          {
            c[i+j*ldc] = alpha * temp;
          }
          else
          {
            c[i+j*ldc] = alpha * temp + beta * c[i+j*ldc];
          }
        }
      }
    }
  }
//
//  Form  C := alpha*A*B' + beta*C
//
  else
  {
    if ( nota )
    {
      for ( j = 0; j < n; j++ )
      {
        if ( beta == 0.0 )
        {
          for ( i = 0; i < m; i++ )
          {
            c[i+j*ldc] = 0.0;
          }
        }
        else if ( beta != 1.0 )
        {
          for ( i = 0; i < m; i++ )
          {
            c[i+j*ldc] = beta * c[i+j*ldc];
          }
        }

        for ( l = 0; l < k; l++ )
        {
          if ( b[j+l*ldb] != 0.0 )
          {
            temp = alpha * b[j+l*ldb];
            for ( i = 0; i < m; i++ )
            {
              c[i+j*ldc] = c[i+j*ldc] + temp * a[i+l*lda];
            }
          }
        }
      }
    }
//
//  Form  C := alpha*A'*B' + beta*C
//
    else
    {
      for ( j = 0; j < n; j++ )
      {
        for ( i = 0; i < m; i++ )
        {
          temp = 0.0;
          for ( l = 0; l < k; l++ )
          {
            temp = temp + a[l+i*lda] * b[j+l*ldb];
          }
          if ( beta == 0.0 )
          {
            c[i+j*ldc] = alpha * temp;
          }
          else
          {
            c[i+j*ldc] = alpha * temp + beta * c[i+j*ldc];
          }
        }
      }
    }
  }

  return;
}
//****************************************************************************80

double franke ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    FRANKE returns the value of the Franke function #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Richard Franke,
//    Scattered Data Interpolation: Tests of Some Methods,
//    Mathematics of Computation,
//    Volume 38, Number 157, January 1982, pages 181-200.
//
//  Parameters:
//
//    Input, double X, Y, the evalution points.
//
//    Output, double FRANKE, the function values.
//
{
  double f;

  f = 
      0.75 * exp ( 
        - ( pow ( 9.0 * x - 2.0, 2 )
          + pow ( 9.0 * y - 2.0, 2 ) ) / 4.0 ) 
    + 0.75 * exp ( 
        - ( pow ( 9.0 * x + 1.0, 2 ) ) / 49.0 
              - ( 9.0 * y + 1.0 )      / 10.0 )       
    + 0.5  * exp ( 
        - ( pow ( 9.0 * x - 7.0, 2 ) 
          + pow ( 9.0 * y - 3.0, 2 ) ) / 4.0 ) 
    - 0.2  * exp ( 
          - pow ( 9.0 * x - 4.0, 2 ) 
          - pow ( 9.0 * y - 7.0, 2 ) );

  return f;
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

void padua2 ( int deg, int degmax, int npd, double wpd[], double fpd[], 
  double raux1[], double raux2[], double c0[], double &esterr )

//****************************************************************************80
//
//  Purpose:
//
//    PADUA2 computes the Padua interpolation coefficient matrix.
//
//  Discussion:
//
//    This function computes the coefficient matrix C0, in the 
//    orthonormal Chebyshev basis T_j(x)T_{k-j}(y), 0 <= j <= k <= DEG, 
//    T_0(x)=1, T_j(x) = sqrt(2) * cos(j * acos(x)), of the 
//    interpolation polynomial of degree DEG of the function values FPD 
//    at the set of NPD Padua points (PD1,PD2) in the square [-1,1]^2. 
//
//    The interpolant may be evaluated at an arbitrary point by the 
//    function PD2VAL. PD1, PD2 and WPD are the Padua points and weights 
//    computed by PDPTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    15 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Parameters:
//
//    Input, int DEG, the degree of approximation.
//
//    Input, int DEGMAX, the maximum degree allowed.
//
//    Input, int NPD, the number of Padua points.
//
//    Input, double WPD[NPD], the weights.
//
//    Input, double FPD[NPD], the value at the Padua points
//    of the function to be interpolated.
//
//    Workspace, double RAUX1[(DEGMAX+1)*(DEG+2)].
//
//    Workspace, double RAUX2[(DEGMAX+1)*(DEG+2)].
//
//    Output, double C0[(DEGMAX+2)*(DEG+1)], the coefficient matrix.
//
//    Output, double &ESTERR, the estimated error.
//
{
  double angle;
  int i;
  int j;
  int k;
  const double pi = 3.1415926535897931;
  double pt;
//
//  Build the matrix P_2 and store it in RAUX2.
//
  for ( i = 0; i <= deg + 1; i++ )
  {
    angle = ( double ) ( i ) * pi / ( double ) ( deg + 1 );
    pt = - cos ( angle );
    cheb ( deg, pt, raux2 + i * ( degmax + 1 ) );
  }
//
//  Build the matrix G(f) and store it in C0.
//
  for ( j = 0; j <= deg + 1; j++ )
  {
    for ( i = 0; i <= degmax + 1; i++ )
    {
      c0[i+j*(degmax+2)] = 0.0;
    }
  }

  k = 0;
  for ( j = 0; j <= deg + 1; j++ )
  {
    for ( i = 0; i <= deg; i++ )
    {
      if ( ( i + j ) % 2 == 0 )
      {
        c0[i+j*(degmax+2)] = fpd[k] * wpd[k];
        k = k + 1;
      }
      else
      {
        c0[i+j*(degmax+2)] = 0.0;
      }
    }
  }
//
//  Compute the matrix-matrix product G(f)*P_2' and store it in RAUX1.
//
  dgemm ( 'n', 't', deg + 1, deg + 1, deg + 2, 1.0, 
    c0, degmax + 2, raux2, degmax + 1, 0.0, raux1, degmax + 1 );
//
//  Build the matrix P_1 and store it in RAUX2.
//
  for ( i = 0; i <= deg; i++ )
  {
    angle = ( double ) ( i ) * pi / ( double ) ( deg );
    pt = - cos ( angle );
    cheb ( deg, pt, raux2 + i * ( degmax + 1 ) );
  }
//
//  Compute the matrix-matrix product C(f) = P_1 * ( G(f) * P_2' ) 
//  and store it in C0.
//
  dgemm ( 'n', 'n', deg + 1, deg + 1, deg + 1, 1.0, 
    raux2, degmax + 1, raux1, degmax + 1, 0.0, c0, degmax + 2 );

  c0[deg+0*(degmax+2)] = c0[deg+0*(degmax+2)] / 2.0;
//
//  Estimate the error.
//
  esterr = 0.0;
  for ( j = 0; j <= 2; j++ )
  {
    for ( i = 0; i <= deg - j; i++ )
    {
      esterr = esterr + fabs ( c0[i+(deg-i-j)*(degmax+2)] );
    }
  }
  esterr = 2.0 * esterr;

  return;
}
//****************************************************************************80

double pd2val ( int deg, int degmax, double c0[], double tg1, double tg2 )

//****************************************************************************80
//
//  Purpose:
//
//    PD2VAL evaluates the Padua2 interpolant.
//
//  Discussion:
//
//    This function returns the value of the interpolant at (TG1,TG2).
//    C0 is the matrix of the coefficients computed by PADUA2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    16 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Parameters:
//
//    Input, int DEG, the degree of approximation.
//
//    Input, int DEGMAX, the maximum degree allowed.         
//
//    Input, double C0[(0:DEGMAX+1)*(0:DEG)], the coefficient matrix.
//
//    Input, double TG1, TG2, the first and second coordinates of
//    the target point.
//
//    Output, double PD2VAL, the value of the interpolant at
//    the target point.
//
{
  int i;
  int j;
  double t;
  double *ttg1;
  double *ttg2;
  double value;
//
//  Compute the normalized Chebyshev polynomials at the target point.
//
  ttg1 = new double[deg+1];
  cheb ( deg, tg1, ttg1 );

  ttg2 = new double[deg+1];
  cheb ( deg, tg2, ttg2 );
//
//  Evaluate the interpolant
//
  value = 0.0;
  for ( i = deg; 0 <= i; i-- )
  {
    t = 0.0;
    for ( j = 0; j <= i; j++ )
    {
      t = t + ttg1[j] * c0[j+(deg-i)*(degmax+2)];
    }
    value = value + ttg2[deg-i] * t;
  }

  delete [] ttg1;
  delete [] ttg2;

  return value;
}
//****************************************************************************80

void pdpts ( int deg, double pd1[], double pd2[], double wpd[], int &npd )

//****************************************************************************80
//
//  Purpose:
//
//    PDPTS returns the points and weights for Padua interpolation.
//
//  Discussion:
//
//    This subroutine computes the first family of Padua points and 
//    weights corresponding to degree DEG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Parameters:
//
//    Input, int DEG, the degree of approximation.
//
//    Output, double PD1[NPD], PD2[NPD], the first and second
//    coordinates of the Padua points
//
//    Output, double WPD[NPD], the weights.
//
//    Output, int &NPD, the number of Padua points.
//    NPD = ( DEG + 1 ) * ( DEG + 2 ) / 2.
//
{
  int itemp0;
  int j;
  int k;
  const double pi = 3.1415926535897931;
  double rtemp0;
//
//  Compute the Padua points of the first family at degree DEG.
//
  if ( deg == 0 )    
  {
    pd1[0] = -1.0;
    pd2[0] = -1.0;
    wpd[0] = 2.0;
    npd = 1;
    return;
  }
   
  npd = 0;
  itemp0 = deg * ( deg + 1 );
  rtemp0 = pi / ( double ) ( itemp0 );

  for ( j = 0; j <= deg + 1; j++ )
  {
    for ( k = ( j % 2 ); k <= deg; k = k + 2 )
    {
      pd1[npd] = - cos ( ( double ) ( ( deg + 1 ) * k ) * rtemp0 );
      pd2[npd] = - cos ( ( double ) ( deg * j ) * rtemp0 );
      wpd[npd] = 2.0 / ( double ) ( itemp0 );

      if ( k == 0 || k == deg )
      {
        wpd[npd] = wpd[npd] / 2.0;
      }

      if ( j == 0 || j == deg + 1 )
      {
        wpd[npd] = wpd[npd] / 2.0;
      }
      npd = npd + 1;
    }
  }

  return;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = - 1.0;
  }
  else
  {
    value = + 1.0;
  }
  return value;
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
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
