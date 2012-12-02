# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <fstream>
# include <cstring>

using namespace std;

# include "box_behnken.hpp"

//****************************************************************************80

double *box_behnken ( int dim_num, int x_num, double range[] )

//****************************************************************************80
//
//  Purpose:
//
//    BOX_BEHNKEN returns a Box-Behnken design for the given number of factors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Box, Donald Behnken,
//    Some new three level designs for the study of quantitative variables,
//    Technometrics,
//    Volume 2, pages 455-475, 1960.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int X_NUM, the number of elements of the design.
//    X_NUM should be equal to DIM_NUM * 2^(DIM_NUM-1) + 1.
//
//    Input, double RANGE[DIM_NUM*2], the minimum and maximum
//    value for each component.
//
//    Output, double BOX_BEHNKEN[DIM_NUM*X_NUM], the elements of the design.
//
{
  int i;
  int i2;
  int j;
  int last_low;
  double *x;
//
//  Ensure that the range is legal.
//
  for ( i = 0; i < dim_num; i++ )
  {
    if ( range[i+1*dim_num] <= range[i+0*dim_num] )
    {
      cerr << "\n";
      cerr << "BOX_BEHNKEN - Fatal error!\n";
      cerr << "  RANGE[" << i << ",1] <= RANGE[" << i << ",0].\n";
      x = NULL;
      exit ( 1 );
    }
  }

  x = new double[dim_num*x_num];
//
//  The first point is the center.
//
  j = 0;

  for ( i = 0; i < dim_num; i++ )
  {
    x[i+j*dim_num] = ( range[i+0*dim_num] + range[i+1*dim_num] ) / 2.0;
  }
//
//  For subsequent elements, one entry is fixed at the middle of the range.
//  The others are set to either extreme.
//
  for ( i = 0; i < dim_num; i++ )
  {
    j = j + 1;
    for ( i2 = 0; i2 < dim_num; i2++ )
    {
      x[i2+j*dim_num] = range[i2+0*dim_num];
    }
    x[i+j*dim_num] = ( range[i+0*dim_num] + range[i+1*dim_num] ) / 2.0;
//
//  The next element is made by finding the last low value, making it
//  high, and all subsequent high values low.
//
    for ( ; ; )
    {
      last_low = -1;

      for ( i2 = 0; i2 < dim_num; i2++ )
      {
        if ( x[i2+j*dim_num] == range[i2+0*dim_num] )
        {
          last_low = i2;
        }
      }

      if ( last_low == -1 )
      {
        break;
      }

      j = j + 1;
      for ( i2 = 0; i2 < dim_num; i2++ )
      {
        x[i2+j*dim_num] = x[i2+(j-1)*dim_num];
      }
      x[last_low+j*dim_num] = range[last_low+1*dim_num];

      for ( i2 = last_low + 1; i2 < dim_num; i2++ )
      {
        if ( x[i2+j*dim_num] == range[i2+1*dim_num] )
        {
          x[i2+j*dim_num] = range[i2+0*dim_num];
        }
      }
    }
  }
  return x;
}
//****************************************************************************80

int box_behnken_size ( int dim_num )

//****************************************************************************80
//
//  Purpose:
//
//    BOX_BEHNKEN_SIZE returns the size of a Box-Behnken design.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Box, Donald Behnken,
//    Some new three level designs for the study of quantitative variables,
//    Technometrics,
//    Volume 2, pages 455-475, 1960.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Output, int X_NUM, the number of elements of the design.
//    X_NUM will be equal to DIM_NUM * 2^(DIM_NUM-1) + 1.
//
{
  int x_num;

  if ( 1 <= dim_num )
  {
    x_num = 1 + dim_num * i4_power ( 2, dim_num - 1 );
  }
  else
  {
    x_num = -1;
  }

  return x_num;
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

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
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

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
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
