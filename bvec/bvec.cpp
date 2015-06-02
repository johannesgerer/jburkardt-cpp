# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "bvec.hpp"

//****************************************************************************80

bool bvec_add ( int n, int bvec1[], int bvec2[], int bvec3[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_ADD adds two (signed) binary vectors.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Example:
//
//    N = 4
//
//      BVEC1       +   BVEC2       =   BVEC3
//
//    ( 0 0 0 0 1 ) + ( 0 0 0 1 1 ) = ( 0 0 1 0 0 )
//
//              1   +           3   =           4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int BVEC1[N], BVEC2[N], the vectors to be added.
//
//    Output, int BVEC3[N], the sum of the two input vectors.
//
//    Output, bool BVEC_ADD, is TRUE if an error occurred.
//
{
  int base = 2;
  int i;
  bool overflow;

  overflow = false;

  for ( i = 0; i < n; i++ )
  {
    bvec3[i] = bvec1[i] + bvec2[i];
  }

  for ( i = n - 1; 0 <= i; i-- )
  {
    while ( base <= bvec3[i] )
    {
      bvec3[i] = bvec3[i] - base;
      if ( 0 < i )
      {
        bvec3[i-1] = bvec3[i-1] + 1;
      }
      else
      {
        overflow = true;
      }
    }
  }

  return overflow;
}
//****************************************************************************80

void bvec_and ( int n, int bvec1[], int bvec2[], int bvec3[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_AND computes the AND of two binary vectors.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int BVEC1[N], BVEC2[N], the vectors.
//
//    Output, int BVEC3[N], the AND of the two vectors.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( bvec1[i] < bvec2[i] )
    {
      bvec3[i] = bvec1[i];
    }
    else
    {
      bvec3[i] = bvec2[i];
    }
  }
  return;
}
//****************************************************************************80

bool bvec_check ( int n, int bvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_CHECK checks a binary vector.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//    The only check made is that the entries are all 0 or 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int BVEC[N], the vector to be checked.
//
//    Output, bool BVEC_CHECK, is TRUE if an error occurred.
//
{
  int base = 2;
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( bvec[i] < 0 || base <= bvec[i] )
    {
      return true;
    }
  }

  return false;
}
//****************************************************************************80

void bvec_complement2 ( int n, int bvec1[], int bvec2[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_COMPLEMENT2 computes the two's complement of a binary vector.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int BVEC1[N], the vector to be two's complemented.
//
//    Output, int BVEC2[N], the two's complemented vector.
//
{
  int base = 2;
  bool error;
  int i;
  int *bvec3;
  int *bvec4;

  bvec3 = new int[n];
  bvec4 = new int[n];
  
  for ( i = 0; i < n; i++ )
  {
    bvec3[i] = ( base - 1 ) - bvec1[i];
  }

  for ( i = 0; i < n - 1; i++ )
  {
    bvec4[i] = 0;
  }
  bvec4[n-1] = 1;

  error = bvec_add ( n, bvec3, bvec4, bvec2 );

  delete [] bvec3;
  delete [] bvec4;
  
  return;
}
//****************************************************************************80

int bvec_enum ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_ENUM enumerates the binary vectors of length N.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.
//
//    BVEC(1) is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
//    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Output, int BVEC_ENUM, the number of binary vectors.
//
{
  int i;
  int value;

  value = 1;
  for ( i = 0; i < n; i++ )
  {
    value = value * 2;
  }

  return value;
}
//****************************************************************************80

void bvec_mul ( int n, int bvec1[], int bvec2[], int bvec3[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_MUL computes the product of two binary vectors.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int BVEC1[N], BVEC2[N], the vectors to be multiplied.
//
//    Output, int BVEC3[N], the product of the two input vectors.
//
{
  int base = 2;
  int *bveca;
  int *bvecb;
  int *bvecc;
  int carry;
  int i;
  int j;
  int product_sign;
//
//  Copy the input.
//
  bveca = new int[n];
  bvecb = new int[n];
  bvecc = new int[n];

  for ( i = 0; i < n; i++ )
  {
    bveca[i] = bvec1[i];
  }

  for ( i = 0; i < n; i++ )
  {
    bvecb[i] = bvec2[i];
  }
//
//  Record the sign of the product.
//  Make the factors positive.
//
  product_sign = 1;

  if ( bveca[0] != 0 )
  {
    product_sign = - product_sign;
    bvec_complement2 ( n, bveca, bveca );
  }

  if ( bvecb[0] != 0 )
  {
    product_sign = - product_sign;
    bvec_complement2 ( n, bvecb, bvecb );
  }

  for ( i = 0; i < n; i++ )
  {
    bvecc[i] = 0;
  }
//
//  Multiply.
//
  for ( i = 1; i <= n - 1; i++ )
  {
    for ( j = 1; j <= n - i; j++ )
    {
      bvecc[j] = bvecc[j] + bveca[n-i] * bvecb[j+i-1];
    }
  }
//
//  Take care of carries.
//
  for ( i = n - 1; 1 <= i; i--)
  {
    carry = bvecc[i] / base;
    bvecc[i] = bvecc[i] - carry * base;
//
//  Unlike the case of BVEC_ADD, we do NOT allow carries into
//  the sign position when multiplying.
//
    if ( 1 < i )
    {
      bvecc[i-1] = bvecc[i-1] + carry;
    }
  }
//
//  Take care of the sign of the product.
//
  if ( product_sign < 0 )
  {
    bvec_complement2 ( n, bvecc, bvecc );
  }
//
//  Copy the output.
//
  for ( i = 0; i < n; i++ )
  {
    bvec3[i] = bvecc[i];
  }

  delete [] bveca;
  delete [] bvecb;
  delete [] bvecc;

  return;
}
//****************************************************************************80

void bvec_next ( int n, int bvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_NEXT generates the next binary vector.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//    The vectors have the order
//
//      (0,0,...,0),
//      (0,0,...,1), 
//      ...
//      (1,1,...,1)
//
//    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
//    we allow wrap around.
//
//  Example:
//
//    N = 3
//
//    Input      Output
//    -----      ------
//    0 0 0  =>  0 0 1
//    0 0 1  =>  0 1 0
//    0 1 0  =>  0 1 1
//    0 1 1  =>  1 0 0
//    1 0 0  =>  1 0 1
//    1 0 1  =>  1 1 0
//    1 1 0  =>  1 1 1
//    1 1 1  =>  0 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input/output, int BVEC[N], on output, the successor to the
//    input vector.  
//
{
  int i;

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( bvec[i] == 0 )
    {
      bvec[i] = 1;
      return;
    }
    bvec[i] = 0;
  }

  return;
}
//****************************************************************************80

void bvec_next_grlex ( int n, int bvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_NEXT generates the next binary vector in GRLEX order.
//
//  Discussion:
//
//    N = 3
//
//    Input      Output
//    -----      ------
//    0 0 0  =>  0 0 1
//    0 0 1  =>  0 1 0
//    0 1 0  =>  1 0 0
//    1 0 0  =>  0 1 1
//    0 1 1  =>  1 0 1
//    1 0 1  =>  1 1 0
//    1 1 0  =>  1 1 1
//    1 1 1  =>  0 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input/output, int BVEC[N], on output, the successor to the
//    input vector.  
//
{
  int i;
  int o;
  int s;
  int z;
//
//  Initialize locations of 0 and 1.
//
  if ( bvec[0] == 0 )
  {
    z = 0;
    o = -1;
  }
  else
  {
    z = -1;
    o = 0;
  }
//
//  Moving from right to left, search for a "1", preceded by a "0".
//
  for ( i = n - 1; 1 <= i; i-- )
  {
    if ( bvec[i] == 1 )
    {
      o = i;
      if ( bvec[i-1] == 0 )
      {
        z = i - 1;
        break;
      }
    }
  }
//
//  BVEC = 0
//
  if ( o == -1 )
  {
    bvec[n-1] = 1;
  }
//
//  01 never occurs.  So for sure, B(1) = 1.
//
  else if ( z == -1 )
  {
    s = 0;
    for ( i = 0; i < n; i++ )
    {
      s = s + bvec[i];
    }
    if ( s == n )
    {
      for ( i = 0; i < n; i++ )
      {
        bvec[i] = 0;
      }
    }
    else
    {
      for ( i = 0; i < n - s - 1; i++ )
      {
        bvec[i] = 0;
      }
      for ( i = n - s - 1; i < n; i++ )
      {
        bvec[i] = 1;
      }
    }
  }
//
//  Found the rightmost "01" string.
//  Replace it by "10".
//  Shift following 1's to the right.
//
  else
  {
    bvec[z] = 1;
    bvec[o] = 0;
    s = 0;
    for ( i = o + 1; i < n; i++ )
    {
      s = s + bvec[i];
    }
    for ( i = o + 1; i < n - s; i++ )
    {
      bvec[i] = 0;
    }
    for ( i = n - s; i < n; i++ )
    {
      bvec[i] = 1;
    }
  }

  return;
}
//****************************************************************************80

void bvec_not ( int n, int bvec1[], int bvec2[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_NOT "negates" or takes the 1's complement of a binary vector.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int BVEC1[N], the vector to be negated.
//
//    Output, int BVEC2[N], the negated vector.
//
{
  int base = 2;
  int i;

  for ( i = 0; i < n; i++ )
  {
    bvec2[i] = ( base - 1 ) - bvec1[i];
  }

  return;
}
//****************************************************************************80

void bvec_or ( int n, int bvec1[], int bvec2[], int bvec3[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_OR computes the OR of two binary vectors.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int BVEC1[N], BVEC2[N], the vectors.
//
//    Output, int BVEC3[N], the OR of the two vectors.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( bvec2[i] < bvec1[i] )
    {
      bvec3[i] = bvec1[i];
    }
    else
    {
      bvec3[i] = bvec2[i];
    }
  }
  return;
}
//****************************************************************************80

void bvec_print ( int n, int bvec[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_PRINT prints a binary integer vector, with an optional title.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int BVEC[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int ihi;
  int ilo;

  if ( 0 < title.length ( ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  for ( ilo = 0; ilo < n; ilo = ilo + 70 )
  {
    ihi = ilo + 70 - 1;
    if ( n - 1 < ihi )
    {
      ihi = n - 1;
    }
    cout << "  ";

    for ( i = ilo; i <= ihi; i++ )
    {
      cout << bvec[i];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void bvec_reverse ( int n, int bvec1[], int bvec2[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_REVERSE reverses a binary vector.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int BVEC1[N], the vector to be reversed.
//
//    Output, int BVEC2[N], the reversed vector.
//
{
  int *bvec3;
  int i;

  bvec3 = new int[n];

  for ( i = 0; i < n; i++ )
  {
    bvec3[i] = bvec1[n-1-i];
  }

  for ( i = 0; i < n; i++ )
  {
    bvec2[i] = bvec3[i];
  }

  free ( bvec3 );

  return;
}
//****************************************************************************80

void bvec_sub ( int n, int bvec1[], int bvec2[], int bvec3[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_SUB subtracts two binary vectors.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Example:
//
//    N = 4
//
//    BVEC1         BVEC2         BVEC3
//    -------       -------       -------
//    0 1 0 0   -   0 0 0 1   =   0 0 1 1
//          4             1             3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int BVEC1[N], BVEC2[N], the vectors to be subtracted.
//
//    Output, int BVEC3[N], the value of BVEC1 - BVEC2.
//
{
  int *bvec4;

  bvec4 = new int[n];

  bvec_complement2 ( n, bvec2, bvec4 );

  bvec_add ( n, bvec1, bvec4, bvec3 );

  delete [] bvec4;

  return;
}
//****************************************************************************80

int bvec_to_i4 ( int n, int bvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_TO_I4 makes an integer from a (signed) binary vector.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Example:
//
//         BVEC   binary  I
//    ----------  -----  --
//    0  0  0  1       1  1
//    0  0  1  0      10  2
//    1  1  0  0    -100 -4
//    0  1  0  0     100  4
//    1  0  0  1    -111 -9
//    1  1  1  1      -0  0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vector.
//
//    Input, int BVEC[N], the binary representation.
//
//    Output, int BVEC_TO_I4, the integer represented by IVEC.
//
{
  int base = 2;
  int *bvec2;
  int i;
  int i4;
  int i4_sign;

  bvec2 = new int[n];

  for ( i = 0; i < n; i++ )
  {
    bvec2[i] = bvec[i];
  }
//
//  Check whether the sign bit is set.
//
  if ( bvec2[0] == base - 1 )
  {
    i4_sign = -1;
    bvec2[0] = 0;
    bvec_complement2 ( n-1, bvec2+1, bvec2+1 );
  }
  else
  {
    i4_sign = 1;
  }

  i4 = 0;
  for ( i = 1; i < n; i++ )
  {
    i4 = base * i4 + bvec2[i];
  }

  i4 = i4_sign * i4;

  delete [] bvec2;

  return i4;
}
//****************************************************************************80

int *bvec_uniform_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_UNIFORM_NEW returns a pseudorandom BVEC.
//
//  Discussion:
//
//    An BVEC is a vector of binary (0/1) values representing an integer.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 December 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre LEcuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the order of the vector.
//
//    Input/output, int &SEED, the "seed" value, which should
//    NOT be 0.  On output, SEED has been updated.
//
//    Output, int BVEC_UNIFORM_NEW[N], a pseudorandom binary vector.
//
{
  int *bvec;
  int i;
  const int i4_huge      = 2147483647;
  const int i4_huge_half = 1073741823;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "BVEC_UNIFORM_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  bvec = new int[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }
    if ( i4_huge_half < seed )
    {
      bvec[i] = 0;
    }
    else
    {
      bvec[i] = 1;
    }
  }
  return bvec;
}
//****************************************************************************80

void bvec_xor ( int n, int bvec1[], int bvec2[], int bvec3[] )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_XOR computes the exclusive OR of two binary vectors.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int BVEC1[N], BVEC2[N], the binary vectors to be XOR'ed.
//
//    Input, int BVEC3[N], the exclusive OR of the two vectors.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    bvec3[i] = ( bvec1[i] + bvec2[i] ) % 2;
  }

  return;
}
//****************************************************************************80

int i4_bclr ( int i4, int pos )

//****************************************************************************80
//
//  Purpose:
//
//    I4_BCLR returns a copy of an I4 in which the POS-th bit is set to 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Military Standard 1753,
//    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
//    9 November 1978.
//
//  Parameters:
//
//    Input, int I4, the integer to be tested.
//
//    Input, int POS, the bit position, between 0 and 31.
//
//    Output, int I4_BCLR, a copy of I4, but with the POS-th bit
//    set to 0.
//
{
  const int i4_huge = 2147483647;
  int j;
  int k;
  int sub;
  int value;

  value = i4;

  if ( pos < 0 )
  {
  }
  else if ( pos < 31 )
  {
    sub = 1;

    if ( 0 <= i4 )
    {
      j = i4;
    }
    else
    {
      j = ( i4_huge + i4 ) + 1;
    }

    for ( k = 1; k <= pos; k++ )
    {
      j = j / 2;
      sub = sub * 2;
    }

    if ( ( j % 2 ) == 1 )
    {
      value = i4 - sub;
    }
  }
  else if ( pos == 31 )
  {
    if ( i4 < 0 )
    {
      value = ( i4_huge + i4 ) + 1;
    }
  }
  else if ( 31 < pos )
  {
    value = i4;
  }

  return value;
}
//****************************************************************************80

int i4_bset ( int i4, int pos )

//****************************************************************************80
//
//  Purpose:
//
//    I4_BSET returns a copy of an I4 in which the POS-th bit is set to 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Military Standard 1753,
//    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
//    9 November 1978.
//
//  Parameters:
//
//    Input, int I4, the integer to be tested.
//
//    Input, int POS, the bit position, between 0 and 31.
//
//    Output, int I4_BSET, a copy of I4, but with the POS-th bit
//    set to 1.
//
{
  int add;
  const int i4_huge = 2147483647;
  int j;
  int k;
  int value;

  value = i4;

  if ( pos < 0 )
  {
  }
  else if ( pos < 31 )
  {
    add = 1;

    if ( 0 <= i4 )
    {
      j = i4;
    }
    else
    {
      j = ( i4_huge + i4 ) + 1;
    }

    for ( k = 1; k <= pos; k++ )
    {
      j = j / 2;
      add = add * 2;
    }

    if ( ( j % 2 ) == 0 )
    {
      value = i4 + add;
    }
  }
  else if ( pos == 31 )
  {
    if ( 0 < i4 )
    {
      value = - ( i4_huge - i4 ) - 1;
    }
  }
  else if ( 31 < pos )
  {
    value = i4;
  }
  return value;
}
//****************************************************************************80

bool i4_btest ( int i4, int pos )

//****************************************************************************80
//
//  Purpose:
//
//    I4_BTEST returns TRUE if the POS-th bit of an I4 is 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Military Standard 1753,
//    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
//    9 November 1978.
//
//  Parameters:
//
//    Input, int I4, the integer to be tested.
//
//    Input, int POS, the bit position, between 0 and 31.
//
//    Output, bool I4_BTEST, is TRUE if the POS-th bit of I4 is 1.
//
{
  const int i4_huge = 2147483647;
  int j;
  int k;
  bool value;

  if ( pos < 0 )
  {
    cerr << "\n";
    cerr << "I4_BTEST - Fatal error!\n";
    cerr << "  POS < 0.\n";
    exit ( 1 );
  }
  else if ( pos < 31 )
  {
    if ( 0 <= i4 )
    {
      j = i4;
    }
    else
    {
      j = ( i4_huge + i4 ) + 1;
    }

    for ( k = 1; k <= pos; k++ )
    {
      j = j / 2;
    }

    if ( ( j % 2 ) == 0 )
    {
      value = false;
    }
    else
    {
      value = true;
    }
  }
  else if ( pos == 31 )
  {
    if ( i4 < 0 )
    {
      value = true;
    }
    else
    {
      value = false;
    }
  }
  else if ( 31 < pos )
  {
    cerr << "\n";
    cerr << "I4_BTEST - Fatal error!\n";
    cerr << "  31 < POS.\n";
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

void i4_to_bvec ( int i4, int n, int bvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_BVEC makes a (signed) binary vector from an I4.
//
//  Discussion:
//
//    A BVEC is a vector of binary digits representing an integer.  
//
//    BVEC[0] is 0 for positive values and 1 for negative values, which
//    are stored in 2's complement form.
//
//    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
//    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//    To guarantee that there will be enough space for any
//    value of I, it would be necessary to set N = 32.
//
//  Example:
//
//    I4       BVEC         binary
//    --  ----------------  ------
//     1  0  0  0  0  0  1      1
//     2  0  0  0  0  1  0     10
//     3  0  0  0  0  1  1     11
//     4  0  0  0  1  0  0    100
//     9  0  0  1  0  0  1   1001
//    -9  1  1  0  1  1  1  -1001 = 110111 (2's complement)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer to be represented.
//
//    Input, int N, the dimension of the vector.
//
//    Output, int BVEC[N], the signed binary representation.
//
{
  int base = 2;
  int i4_copy;
  int j;

  i4_copy = abs ( i4 );

  for ( j = n - 1; 1 <= j; j-- )
  {
    bvec[j] = ( i4_copy % base );

    i4_copy = i4_copy / base;
  }

  bvec[0] = 0;

  if ( i4 < 0 )
  {
    bvec_complement2 ( n, bvec, bvec );
  }

  return;
}
//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
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
