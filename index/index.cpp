# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "index.hpp"

//****************************************************************************80

int index0 ( int i_min, int i, int i_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX0 indexes a 1D vector using a zero base.
//
//  Discussion:
//
//    Index       Element
//    ---------   --------
//    0           I_MIN
//    INDEX0      I
//   (INDEX_MAX)  I_MAX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for the first index,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX0, the index of element I.
//
{
  int index_min = 0;
  int value;

  value = index_min + ( i - i_min );

  return value;
}
//****************************************************************************80

int index01 ( int i_min, int i, int i_max, int j_min, int j, int j_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX01 indexes a 2D array by columns, with a zero base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
//    and increasing the row index first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX01, the index of element (I,J).
//
{
  int index_min = 0;
  int value;

  value = 
    index_min 
    + (         i - i_min ) 
    + ( i_max + 1 - i_min ) * ( j - j_min );

  return value;
}
//****************************************************************************80

int index012 ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int k_min, int k, int k_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX012 indexes a 3D array by columns with zero base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry (I_MIN,J_MIN,K_MIN), 
//    and increasing the row index first, then the column index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Input, int K_MIN, K, K_MAX, for plane indices,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX012, the index of element (I,J,K).
//
{
  int index_min = 0;
  int value;

  value = 
      index_min 
    + (         i - i_min ) 
    + ( i_max + 1 - i_min ) * (         j - j_min ) *  
    + ( i_max + 1 - i_min ) * ( j_max + 1 - j_min ) * ( k - k_min );

  return value;
}
//****************************************************************************80

int index0123 ( int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max, 
  int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX0123 indexes a 4D array by columns, with a zero base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
//    and increasing the initial index first, then the second, third and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1_MIN, I1, I1_MAX, for index 1,
//    the minimum, the index, and the maximum.
//
//    Input, int I2_MIN, I2, I2_MAX, for index 2,
//    the minimum, the index, and the maximum.
//
//    Input, int I3_MIN, I3, I3_MAX, for index 3,
//    the minimum, the index, and the maximum.
//
//    Input, int I4_MIN, I4, I4_MAX, for index 4,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX0123, the index of (I1,I2,I3,I4).
//
{
  int index_min = 0;
  int value;

  value = 
      index_min 
    + (         i1 - i1_min ) 
    + ( i1_max + 1 - i1_min ) * (         i2 - i2_min ) 
    + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) 
    * (         i3 - i3_min ) 
    + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) 
    * ( i3_max + 1 - i3_min ) * (         i4 - i4_min );

  return value;
}
//****************************************************************************80

int index0n ( int n, int i_min[], int i[], int i_max[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX0N indexes an N-dimensional array by columns, with zero base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry 
//      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
//    and increasing the first index up to I_MAX(1), 
//    then the second and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of indices.
//
//    Input, int I_MIN[N], the minimum indices.
//
//    Input, int I[N], the indices.
//
//    Input, int I_MAX[N], for maximum indices.
//
//    Output, int INDEX0N, the index of element I.
//
{
  int index_min = 0;
  int j;
  int value;

  value = ( i[n-1] - i_min[n-1] );

  for ( j = n - 2; 0 <= j; j-- )
  {
    value = value * ( i_max[j] + 1 - i_min[j] ) + ( i[j] - i_min[j] );
  }
  value = value + index_min;

  return value;
}
//****************************************************************************80

int index1 ( int i_min, int i, int i_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX1 indexes a 1D vector using a unit base.
//
//  Discussion:
//
//    Index       Element
//    ---------   --------
//    1           I_MIN
//    INDEX1      I
//   (INDEX_MAX)  I_MAX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for the first index,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX1, the index of element I.
//
{
  int index_min = 1;
  int value;

  value = index_min + ( i - i_min );

  return value;
}
//****************************************************************************80

int index10 ( int i_min, int i, int i_max, int j_min, int j, int j_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX10 indexes a 2D array by rows, with a zero base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
//    and increasing the column index first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX10, the index of element (I,J).
//
{
  int index_min = 0;
  int value;

  value = index_min 
             +                         ( j - j_min ) 
             + ( i - i_min ) * ( j_max + 1 - j_min );

  return value;
}
//****************************************************************************80

int index12 ( int i_min, int i, int i_max, int j_min, int j, int j_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX12 indexes a 2D array by columns, with a unit base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
//    and increasing the row index first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX12, the index of element (I,J).
//
{
  int index_min = 1;
  int value;

  value = 
    index_min 
    + (         i - i_min ) 
    + ( i_max + 1 - i_min ) * ( j - j_min );

  return value;
}
//****************************************************************************80

int index123 ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int k_min, int k, int k_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX123 indexes a 3D array by columns with unit base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry (I_MIN,J_MIN,K_MIN), 
//    and increasing the row index first, then the column index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Input, int K_MIN, K, K_MAX, for plane indices,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX123, the index of element (I,J,K).
//
{
  int index_min = 1;
  int value;

  value = 
      index_min 
    + (         i - i_min ) 
    + ( i_max + 1 - i_min ) * (         j - j_min ) *  
    + ( i_max + 1 - i_min ) * ( j_max + 1 - j_min ) * ( k - k_min );

  return value;
}
//****************************************************************************80

int index1234 ( int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max, 
  int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX1234 indexes a 4D array by columns, with a unit base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
//    and increasing the initial index first, then the second, third and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1_MIN, I1, I1_MAX, for index 1,
//    the minimum, the index, and the maximum.
//
//    Input, int I2_MIN, I2, I2_MAX, for index 2,
//    the minimum, the index, and the maximum.
//
//    Input, int I3_MIN, I3, I3_MAX, for index 3,
//    the minimum, the index, and the maximum.
//
//    Input, int I4_MIN, I4, I4_MAX, for index 4,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX1234, the index of (I1,I2,I3,I4).
//
{
  int index_min = 1;
  int value;

  value = 
      index_min 
    + (         i1 - i1_min ) 
    + ( i1_max + 1 - i1_min ) * (         i2 - i2_min ) 
    + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) 
    * (         i3 - i3_min ) 
    + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) 
    * ( i3_max + 1 - i3_min ) * (         i4 - i4_min );

  return value;
}
//****************************************************************************80

int index1n ( int n, int i_min[], int i[], int i_max[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX1N indexes an N-dimensional array by columns, with unit base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry 
//      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
//    and increasing the first index up to I_MAX(1), 
//    then the second and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of indices.
//
//    Input, int I_MIN[N], the minimum indices.
//
//    Input, int I[N], the indices.
//
//    Input, int I_MAX[N], for maximum indices.
//
//    Output, int INDEX1N, the index of element I.
//
{
  int index_min = 1;
  int j;
  int value;

  value = ( i[n-1] - i_min[n-1] );

  for ( j = n - 2; 0 <= j; j-- )
  {
    value = value * ( i_max[j] + 1 - i_min[j] ) + ( i[j] - i_min[j] );
  }
  value = value + index_min;

  return value;
}
//****************************************************************************80

int index21 ( int i_min, int i, int i_max, int j_min, int j, int j_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX21 indexes a 2D array by rows, with a unit base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
//    and increasing the column index first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX21, the index of element (I,J).
//
{
  int index_min = 1;
  int value;

  value = index_min 
             +                         ( j - j_min ) 
             + ( i - i_min ) * ( j_max + 1 - j_min );

  return value;
}
//****************************************************************************80

int index210 ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int k_min, int k, int k_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX210 indexes a 3D array by rows, with zero base.
//
//  Discussion:
//
//    When we say "by rows", we really just mean that entries of the array are 
//    indexed starting at entry (I_MIN,J_MIN,K_MIN), and the increasing the LAST
//    index first, then the next-to-the-last, and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Input, int K_MIN, K, K_MAX, for plane indices,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX210, the index of element (I,J,K).
//
{
  int index_min = 0;
  int value;

  value = 
      index_min 
    +                                                 ( k - k_min ) 
    +                         ( j - j_min ) * ( k_max + 1 - k_min ) 
    + ( i - i_min ) * ( j_max + 1 - j_min ) * ( k_max + 1 - k_min );

  return value;
}
//****************************************************************************80

int index321 ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int k_min, int k, int k_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX321 indexes a 3D array by rows, with zero base.
//
//  Discussion:
//
//    When we say "by rows", we really just mean that entries of the array are 
//    indexed starting at entry (I_MIN,J_MIN,K_MIN), and the increasing the LAST
//    index first, then the next-to-the-last, and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Input, int K_MIN, K, K_MAX, for plane indices,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX321, the index of element (I,J,K).
//
{
  int index_min = 1;
  int value;

  value = 
      index_min 
    +                                                 ( k - k_min ) 
    +                         ( j - j_min ) * ( k_max + 1 - k_min ) 
    + ( i - i_min ) * ( j_max + 1 - j_min ) * ( k_max + 1 - k_min );

  return value;
}
//****************************************************************************80

int index3210 ( int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max, 
  int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX3210 indexes a 4D array by rows, with zero base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
//    and increasing the last index, then the next to last, and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1_MIN, I1, I1_MAX, for index 1,
//    the minimum, the index, and the maximum.
//
//    Input, int I2_MIN, I2, I2_MAX, for index 2,
//    the minimum, the index, and the maximum.
//
//    Input, int I3_MIN, I3, I3_MAX, for index 3,
//    the minimum, the index, and the maximum.
//
//    Input, int I4_MIN, I4, I4_MAX, for index 4,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX3210, the index of (I1,I2,I3,I4).
//
{
  int index_min = 0;
  int value;

  value = 
    index_min 
    +         ( i4 - i4_min ) 
    +                                                     ( i3 - i3_min ) 
    * ( i4_max + 1 - i4_min ) 
    +                           ( i2 - i2_min ) * ( i3_max + 1 - i3_min ) 
    * ( i4_max + 1 - i4_min ) 
    + ( i1 - i1_min ) * ( i2_max + 1 - i2_min ) * ( i3_max + 1 - i3_min ) 
    * ( i4_max + 1 - i4_min );

  return value;
}
//****************************************************************************80

int index4321 ( int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max, 
  int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX4321 indexes a 4D array by rows, with unit base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
//    and increasing the last index, then the next to last, and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1_MIN, I1, I1_MAX, for index 1,
//    the minimum, the index, and the maximum.
//
//    Input, int I2_MIN, I2, I2_MAX, for index 2,
//    the minimum, the index, and the maximum.
//
//    Input, int I3_MIN, I3, I3_MAX, for index 3,
//    the minimum, the index, and the maximum.
//
//    Input, int I4_MIN, I4, I4_MAX, for index 4,
//    the minimum, the index, and the maximum.
//
//    Output, int INDEX4321, the index of (I1,I2,I3,I4).
//
{
  int index_min = 1;
  int value;

  value = 
    index_min 
    +         ( i4 - i4_min ) 
    +                                                     ( i3 - i3_min ) 
    * ( i4_max + 1 - i4_min ) 
    +                           ( i2 - i2_min ) * ( i3_max + 1 - i3_min ) 
    * ( i4_max + 1 - i4_min ) 
    + ( i1 - i1_min ) * ( i2_max + 1 - i2_min ) * ( i3_max + 1 - i3_min ) 
    * ( i4_max + 1 - i4_min );

  return value;
}
//****************************************************************************80

int indexn0 ( int n, int i_min[], int i[], int i_max[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEXN0 indexes an N-dimensional array by rows, with zero base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry 
//      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
//    and increasing the last index up to I_MAX(N), 
//    then the next-to-last and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of indices.
//
//    Input, int I_MIN[N], the minimum indices.
//
//    Input, int I[N], the indices.
//
//    Input, int I_MAX[N], for maximum indices.
//
//    Output, int INDEXN0, the index of element I.
//
{
  int index_min = 0;
  int j;
  int value;

  value = ( i[0] - i_min[0] );

  for ( j = 1; j < n; j++ )
  {
    value = value * ( i_max[j] + 1 - i_min[j] ) + ( i[j] - i_min[j] );
  }
  value = value + index_min;

  return value;
}
//****************************************************************************80

int indexn1 ( int n, int i_min[], int i[], int i_max[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEXN1 indexes an N-dimensional array by rows, with unit base.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry 
//      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
//    and increasing the last index up to I_MAX(N), 
//    then the next-to-last and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of indices.
//
//    Input, int I_MIN[N], the minimum indices.
//
//    Input, int I[N], the indices.
//
//    Input, int I_MAX[N], for maximum indices.
//
//    Output, int INDEXN1, the index of element I.
//
{
  int index_min = 1;
  int j;
  int value;

  value = ( i[0] - i_min[0] );

  for ( j = 1; j < n; j++ )
  {
    value = value * ( i_max[j] + 1 - i_min[j] ) + ( i[j] - i_min[j] );
  }
  value = value + index_min;

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

