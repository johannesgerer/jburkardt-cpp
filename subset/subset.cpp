# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "subset.hpp"

//****************************************************************************80

int asm_enum ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    ASM_ENUM returns the number of alternating sign matrices of a given order.
//
//  Discussion:
//
//    N     ASM_NUM
//
//    0       1
//    1       1
//    2       2
//    3       7
//    4      42
//    5     429
//    6    7436
//    7  218348
//
//    A direct formula is
//
//      ASM_NUM ( N ) = product ( 0 <= I <= N-1 ) ( 3 * I + 1 )! / ( N + I )!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrices.
//
//    Output, int ASM_ENUM, the number of alternating sign matrices of
//    order N.
//
{
  int *a;
  int asm_num;
  int *b;
  int *c;
  int i;
  int nn;

  if ( n + 1 <= 0 )
  {
    return 0;
  }
//
//  Row 1
//
  if ( n + 1 == 1 )
  {
    return 1;
  }
//
//  Row 2
//
  if ( n + 1 == 2 )
  {
    return 1;
  }

  a = new int[n+1];
  b = new int[n+1];
  c = new int[n+1];

  b[0] = 2;
  c[0] = 2;
  a[0] = 1;
  a[1] = 1;
//
//  Row 3 and on.
//
  for ( nn = 3; nn <= n; nn++ )
  {
    b[nn-2] = nn;
    for ( i = nn-2; 2 <= i; i-- )
    {
      b[i-1] = b[i-1] + b[i-2];
    }
    b[0] = 2;

    c[nn-2] = 2;
    for ( i = nn-2; 2 <= i; i--)
    {
      c[i-1] = c[i-1] + c[i-2];
    }
    c[0] = nn;

    for ( i = 2; i <= nn-1; i++ )
    {
      a[0] = a[0] + a[i-1];
    }

    for  ( i = 2; i <= nn; i++ )
    {
      a[i-1] = a[i-2] * c[i-2] / b[i-2];
    }
  }

  asm_num = 0;
  for ( i = 0; i < n; i++ )
  {
    asm_num = asm_num + a[i];
  }

  delete [] a;
  delete [] b;
  delete [] c;

  return asm_num;
}
//****************************************************************************80

void asm_triangle ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    ASM_TRIANGLE returns a row of the alternating sign matrix triangle.
//
//  Discussion:
//
//    The first seven rows of the triangle are as follows:
//
//          1      2      3      4      5      6     7
//
//    0     1
//    1     1      1
//    2     2      3      2
//    3     7     14     14      7
//    4    42    105    135    105     42
//    5   429   1287   2002   2002   1287    429
//    6  7436  26026  47320  56784  47320  26026  7436
//
//    For a given N, the value of A(J) represents entry A(I,J) of
//    the triangular matrix, and gives the number of alternating sign matrices
//    of order N in which the (unique) 1 in row 1 occurs in column J.
//
//    Thus, of alternating sign matrices of order 3, there are
//    2 with a leading 1 in column 1:
//
//      1 0 0  1 0 0
//      0 1 0  0 0 1
//      0 0 1  0 1 0
//
//    3 with a leading 1 in column 2, and
//
//      0 1 0  0 1 0  0 1 0
//      1 0 0  0 0 1  1-1 1
//      0 0 1  1 0 0  0 1 0
//
//    2 with a leading 1 in column 3:
//
//      0 0 1  0 0 1
//      1 0 0  0 1 0
//      0 1 0  1 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the desired row.
//
//    Output, int A[N+1], the entries of the row.
//
{
  int *b;
  int *c;
  int i;
  int nn;
//
  if ( n + 1 <= 0 )
  {
    return;
  }
//
//  Row 1
//
  a[0] = 1;

  if ( n + 1 == 1 )
  {
    return;
  }
//
//  Row 2
//
  a[0] = 1;
  a[1] = 1;

  if ( n + 1 == 2 )
  {
    return;
  }
//
//  Row 3 and on.
//
  b = new int[n+1];
  c = new int[n+1];

  b[0] = 2;
  c[0] = 2;

  for ( nn = 3; nn <= n+1; nn++ )
  {

    b[nn-2] = nn;
    for ( i = nn-2; 2 <= i; i-- )
    {
      b[i-1] = b[i-1] + b[i-2];
    }
    b[0] = 2;

    c[nn-2] = 2;
    for ( i = nn-2; 2 <= i; i-- )
    {
      c[i-1] = c[i-1] + c[i-2];
    }
    c[0] = nn;

    for ( i = 2; i <= nn-1; i++ )
    {
      a[0] = a[0] + a[i-1];
    }

    for ( i = 2; i <= nn; i++ )
    {
      a[i-1] = a[i-2] * c[i-2] / b[i-2];
    }

  }

  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void bell ( int n, int b[] )

//****************************************************************************80
//
//  Purpose:
//
//    BELL returns the Bell numbers from 0 to N.
//
//  Discussion:
//
//    The Bell number B(N) is the number of restricted growth functions
//    on N.
//
//    Note that the Stirling numbers of the second kind, S^m_n, count the
//    number of partitions of N objects into M classes, and so it is
//    true that
//
//      B(N) = S^1_N + S^2_N + ... + S^N_N.
//
//    The Bell number B(N) is defined as the number of partitions (of
//    any size) of a set of N distinguishable objects.
//
//    A partition of a set is a division of the objects of the set into
//    subsets.
//
//    For instance, there are 15 partitions of a set of 4 objects:
//
//      (1234), (123)(4), (124)(3), (12)(34), (12)(3)(4),
//      (134)(2), (13)(24), (13)(2)(4), (14)(23), (1)(234),
//      (1)(23)(4), (14)(2)(3), (1)(24)(3), (1)(2)(34), (1)(2)(3)(4)
//
//    and so B(4) = 15.
//
//    The recursion formula is:
//
//      B(I) = sum ( 1 <= J <= I ) Binomial ( I-1, J-1 ) * B(I-J)
//
//  Example:
//
//     N         B(N)
//     0           1
//     1           1
//     2           2
//     3           5
//     4          15
//     5          52
//     6         203
//     7         877
//     8        4140
//     9       21147
//    10      115975
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of Bell numbers desired.
//
//    Output, int B[N+1], the Bell numbers from 0 to N.
//
{
  int i;
  int j;

  b[0] = 1;

  for ( i = 1; i <= n; i++ )
  {
    b[i] = 0;
    for ( j = 1; j <= i; j++ )
    {
      b[i] = b[i] + b[i-j] * i4_choose ( i-1, j-1 );
    }
  }

  return;
}
//****************************************************************************80

void bell_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    BELL_VALUES returns some values of the Bell numbers for testing.
//
//  Modified:
//
//    08 May 2003
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int &N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int &N, the order of the Bell number.
//
//    Output, int &C, the value of the Bell number.
//
{
# define N_MAX 11

  int c_vec[N_MAX] = { 1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 };
  int n_vec[N_MAX] = { 0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10};

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  if ( N_MAX <= n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data];
    c = c_vec[n_data];
    n_data = n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void binary_vector_next ( int n, int bvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    BINARY_VECTOR_NEXT generates the next binary vector.
//
//  Discussion:
//
//    A binary vector is a vector whose entries are 0 or 1.
//
//    The user inputs an initial zero vector to start.  The program returns
//    the "next" vector.
//
//    The vectors are produced in the order:
//
//    ( 0, 0, 0, ..., 0 )
//    ( 1, 0, 0, ..., 0 ) 
//    ( 0, 1, 0, ..., 0 )
//    ( 1, 1, 0, ..., 0 )
//    ( 0, 0, 1, ..., 0 )
//    ( 1, 0, 1, ..., 0 )
//               ...
//    ( 1, 1, 1, ..., 1)
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
//    0 0 0  =>  1 0 0
//    1 0 0  =>  0 1 0
//    0 1 0  =>  1 1 0
//    1 1 0  =>  0 0 1
//    0 0 1  =>  1 0 1
//    1 0 1  =>  0 1 1
//    0 1 1  =>  1 1 1
//    1 1 1  =>  0 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input/output, int BVEC[N], on output, the successor 
//    to the input vector.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {  
    if ( bvec[i] == 1 )
    {
      bvec[i] = 0;
    }
    else 
    {
      bvec[i] = 1;
      break;
    }
  }
  return;
}
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
    bvec3[i] = i4_min ( bvec1[i], bvec2[i] );
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
    bvec3[i] = i4_max ( bvec1[i], bvec2[i] );
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

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  for ( ilo = 0; ilo < n; ilo = ilo + 70 )
  {
    ihi = i4_min ( ilo + 70 - 1, n - 1 );
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

  for ( i = 0; i < n; i++ )
  {
    bvec3[i] = bvec1[n-1-i];
  }

  for ( i = 0; i < n; i++ )
  {
    bvec2[i] = bvec3[i];
  }

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
  int i;

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

void catalan ( int n, int c[] )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN computes the Catalan numbers, from C(0) to C(N).
//
//  Discussion:
//
//    The Catalan number C(N) counts:
//
//    1) the number of binary trees on N vertices;
//    2) the number of ordered trees on N+1 vertices;
//    3) the number of full binary trees on 2N+1 vertices;
//    4) the number of well formed sequences of 2N parentheses;
//    5) the number of ways 2N ballots can be counted, in order,
//       with N positive and N negative, so that the running sum
//       is never negative;
//    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
//    7) the number of monotone functions from [1..N} to [1..N} which
//       satisfy f(i) <= i for all i;
//    8) the number of ways to triangulate a polygon with N+2 vertices.
//
//    The formula is:
//
//      C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
//           = 1 / (N+1) * COMB ( 2N, N )
//           = 1 / (2N+1) * COMB ( 2N+1, N+1).
//
//  First values:
//
//     C(0)     1
//     C(1)     1
//     C(2)     2
//     C(3)     5
//     C(4)    14
//     C(5)    42
//     C(6)   132
//     C(7)   429
//     C(8)  1430
//     C(9)  4862
//    C(10) 16796
//
//  Recursion:
//
//    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
//    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
//
//  Example:
//
//    N = 3
//
//    ()()()
//    ()(())
//    (()())
//    (())()
//    ((()))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input, int N, the number of Catalan numbers desired.
//
//    Output, int C[N+1], the Catalan numbers from C(0) to C(N).
//
{
  int i;

  if ( n < 0 )
  {
    return;
  }

  c[0] = 1;
//
//  The extra parentheses ensure that the integer division is
//  done AFTER the integer multiplication.
//
  for ( i = 1; i <= n; i++ )
  {
    c[i] = ( c[i-1] * 2 * ( 2 * i - 1 ) ) / ( i + 1 );
  }

  return;
}
//****************************************************************************80

void catalan_row_next ( bool next, int n, int irow[] )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN_ROW computes row N of Catalan's triangle.
//
//  Example:
//
//    I\J 0   1   2   3   4   5   6
//
//    0   1
//    1   1   1
//    2   1   2   2
//    3   1   3   5   5
//    4   1   4   9  14  14
//    5   1   5  14  28  42  42
//    6   1   6  20  48  90 132 132
//
//  Recursion:
//
//    C(0,0) = 1
//    C(I,0) = 1
//    C(I,J) = 0 for I < J
//    C(I,J) = C(I,J-1) + C(I-1,J)
//    C(I,I) is the I-th Catalan number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, bool NEXT, indicates whether this is a call for
//    the 'next' row of the triangle.
//    NEXT = FALSE, this is a startup call.  Row N is desired, but
//    presumably this is a first call, or row N-1 was not computed
//    on the previous call.
//    NEXT = TRUE, this is not the first call, and row N-1 was computed
//    on the previous call.  In this case, much work can be saved
//    by using the information from the previous values of IROW
//    to build the next values.
//
//    Input, int N, the index of the row of the triangle desired.
//
//    Input/output, int IROW[N+1], the row of coefficients.
//    If NEXT = FALSE, then IROW is not required to be set on input.
//    If NEXT = TRUE, then IROW must be set on input to the value of
//    row N-1.
//
{
  int i;
  int j;
//
  if ( n < 0 )
  {
    return;
  }

  if ( !next )
  {
    irow[0] = 1;
    for ( i = 1; i <= n; i++ )
    {
      irow[i] = 0;
    }

    for ( i = 1; i <= n; i++ )
    {
      irow[0] = 1;

      for ( j = 1; j <= i-1; j++ )
      {
        irow[j] = irow[j] + irow[j-1];
      }

      irow[i] = irow[i-1];

    }
  }
  else
  {
    irow[0] = 1;

    for ( j = 1; j <= n-1; j++ )
    {
      irow[j] = irow[j] + irow[j-1];
    }

    if ( 1 <= n )
   {
      irow[n] = irow[n-1];
    }

  }

  return;
}
//****************************************************************************80

void catalan_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN_VALUES returns some values of the Catalan numbers for testing.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int &N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int &N, the order of the Catalan number.
//
//    Output, int &C, the value of the Catalan number.
//
{
# define N_MAX 11

  int c_vec[N_MAX] = { 1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796 };
  int n_vec[N_MAX] = { 0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  if ( N_MAX <= n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data];
    c = c_vec[n_data];
    n_data = n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cbt_traverse ( int depth )

//****************************************************************************80
//
//  Purpose:
//
//    CBT_TRAVERSE traverses a complete binary tree of given depth.
//
//  Discussion:
//
//    There will be 2^DEPTH terminal nodes of the complete binary tree.
//
//    This function traverses the tree, and prints out a binary code of 0's
//    and 1's each time it encounters a terminal node.  This results in a 
//    printout of the binary digits from 0 to 2^DEPTH - 1.
//
//    The function is intended as a framework to be used to traverse a binary
//    tree.  Thus, in practice, a user would insert some action when a terminal
//    node is encountered.
//
//    Another use would occur when a combinatorial search is being made, for
//    example in a knapsack problem.  Each binary string then represents which
//    objects are to be included in the knapsack.  In that case, the traversal
//    could be speeded up by noticing cases where a nonterminal node has been
//    reached, but the knapsack is already full, in which case the only solution
//    uses none of the succeeding items, or overfull, in which case no solutions
//    exist that include this initial path segment.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DEPTH, the depth of the tree.
//
{
  int *b;
  int direction;
  int DOWNLEFT = 1;
  int i;
  int k;
  int p;
  int UP = 3;
  int UPDOWNRIGHT = 2;

  if ( depth < 1 )
  {
    return;
  }

  b = new int[depth];

  for ( i = 0; i < depth; i++ )
  {
    b[i] = 0;
  }
  p = 0;
  direction = DOWNLEFT;
  k = 0;

  for ( ; ; )
  {
//
//  Try going in direction DOWNLEFT.
//
    if ( direction == DOWNLEFT )
    {
      p = p + 1;
      b[p] = 0;
      if ( depth <= p )
      {
        cout << "  " << setw(4) << k;
        for ( i = 0; i < depth; i++ )
        {
          cout << setw(1) << b[i];
        }
        cout << "\n";
        k = k + 1;
        direction = UPDOWNRIGHT;
      }
    }
//
//  Try going in direction UPDOWNRIGHT.
//
    if ( direction == UPDOWNRIGHT )
    {
      b[p] = + 1;
      if ( p < depth )
      {
        direction = DOWNLEFT;
      }
      else
      {
        cout << "  " << setw(4) << k;
        for ( i = 0; i < depth; i++ )
        {
          cout << setw(1) << b[i];
        }
        cout << "\n";
        k = k + 1;
        direction = UP;
      }
    }
//
//  Try going in direction UP.
//
    if ( direction == UP )
    {
      p = p - 1;
      if ( 1 <= p )
      {
        if ( b[p] == 0 )
        {
          direction = UPDOWNRIGHT;
        }
      }
      else
      {
        break;
      }
    }
  }

  delete [] b;

  return;
}
//****************************************************************************80

void cfrac_to_rat ( int n, int a[], int p[], int q[] )

//****************************************************************************80
//
//  Purpose:
//
//    CFRAC_TO_RAT converts a monic continued fraction to an ordinary fraction.
//
//  Discussion:
//
//    The routine is given the monic or "simple" continued fraction with
//    integer coefficients:
//
//      A(1) + 1 / ( A(2) + 1 / ( A(3) ... + 1 / A(N) ) )
//
//    and returns the N successive approximants P(I)/Q(I)
//    to the value of the rational number represented by the continued
//    fraction, with the value exactly equal to the final ratio P(N)/Q(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2004
//
//  Author:
//
//    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
//    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, 
//    Christoph Witzgall.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
//    John Rice, Henry Thatcher, Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968.
//
//  Parameters:
//
//    Input, int N, the number of continued fraction coefficients.
//
//    Input, int A[N], the continued fraction coefficients.
//
//    Output, int P[N], Q[N], the N successive approximations
//    to the value of the continued fraction.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      p[i] = a[i] * 1 + 0;
      q[i] = a[i] * 0 + 1;
    }
    else if ( i == 1 )
    {
      p[i] = a[i] * p[i-1] + 1;
      q[i] = a[i] * q[i-1] + 0;
    }
    else
    {
      p[i] = a[i] * p[i-1] + p[i-2];
      q[i] = a[i] * q[i-1] + q[i-2];
    }

  }

  return;
}
//****************************************************************************80

void cfrac_to_rfrac ( int m, double g[], double h[], double p[], double q[] )

//****************************************************************************80
//
//  Purpose:
//
//    CFRAC_TO_RFRAC converts a polynomial fraction from continued to rational form.
//
//  Discussion:
//
//    The routine accepts a continued polynomial fraction:
//
//      G(1)     / ( H(1) +
//      G(2) * X / ( H(2) +
//      G(3) * X / ( H(3) + ...
//      G(M) * X / ( H(M) )...) ) )
//
//    and returns the equivalent rational polynomial fraction:
//
//      P(1) + P(2) * X + ... + P(L1) * X**(L1)
//      -------------------------------------------------------
//      Q(1) + Q(2) * X + ... + Q(L2) * X**(L2-1)
//
//    where
//
//      L1 = (M+1)/2
//      L2 = (M+2)/2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2004
//
//  Author:
//
//    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
//    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, 
//    Christoph Witzgall.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
//    John Rice, Henry Thatcher, Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968.
//
//  Parameters:
//
//    Input, int M, the number of continued fraction polynomial coefficients.
//
//    Input, double G[M], H[M], the continued polynomial fraction coefficients.
//
//    Output, double P[(M+1)/2], Q[(M+2)/2], the rational polynomial fraction
//    coefficients.
//
{
  double *a;
  int i;
  int j;

  if ( m == 1 )
  {
    p[0] = g[0];
    q[0] = h[0];
    return;
  }

  a = new double[m*((m+2)/2)];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < (m+2)/2; j++ )
    {
      a[i+j*m] = 0.0;
    }
  }
//
//  Solve for P's
//
  a[0+0*m] = g[0];
  a[1+0*m] = g[0] * h[1];

  for ( i = 2; i < m; i++ )
  {
    a[i+0*m] = h[i] * a[i-1+0*m];
    for ( j = 1; j < (i+2)/2; j++ )
    {
      a[i+j*m] = h[i] * a[i-1+j*m] + g[i] * a[i-2+(j-1)*m];
    }
  }

  for ( j = 0; j < (m+1)/2; j++ )
  {
    p[j] = a[m-1+j*m];
  }
//
//  Solve for Q's.
//
  a[0+0*m] = h[0];
  a[1+0*m] = h[0] * h[1];
  a[1+1*m] = g[1];

  for ( i = 2; i < m; i++ )
  {
    a[i+0*m] = h[i] * a[i-1+0*m];
    for ( j = 1; j < (i+3)/2; j++ )
    {
      a[i+j*m] = h[i] * a[i-1+j*m] + g[i] * a[i-2+(j-1)*m];
    }
  }

  for ( j = 0; j < (m+2)/2; j++ )
  {
    q[j] = a[m-1+j*m];
  }

  delete [] a;

  return;
}
//****************************************************************************80

char ch_cap ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= c && c <= 122 ) 
  {
    c = c - 32;
  }   

  return c;
}
//****************************************************************************80

void change_greedy ( int total, int coin_num, int coin_value[], int &change_num, 
  int change[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHANGE_GREEDY makes change for a given total using the biggest coins first.
//
//  Discussion:
//
//    The algorithm is simply to use as many of the largest coin first,
//    then the next largest, and so on.
//
//    It is assumed that there is always a coin of value 1.  The
//    algorithm will otherwise fail!
//
//  Example:
//
//    Total = 17
//    COIN_NUM = 3
//    COIN_VALUE = (/ 1, 5, 10 /)
//
//
//    #  CHANGE              COIN_VALUE(CHANGE)
//
//    4  3 2 1 1             10 5 1 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TOTAL, the total for which change is to be made.
//
//    Input, int COIN_NUM, the number of types of coins.
//
//    Input, int COIN_VALUE[COIN_NUM], the value of each coin.
//    The values should be in ascending order, and if they are not,
//    they will be sorted.
//
//    Output, int &CHANGE_NUM, the number of coins given in change.
//
//    Output, int CHANGE[TOTAL], the indices of the coins will be
//    in entries 1 through CHANGE_NUM.
//
{
  int j;

  change_num = 0;
//
//  Find the largest coin smaller than the total.
//
  j = coin_num - 1;

  while ( 0 <= j )
  {
    if ( coin_value[j] <= total )
    {
      break;
    }
    j = j - 1;
  }

  if ( j < 0 )
  {
    return;
  }
//
//  Subtract the current coin from the total.
//  Once that coin is too big, use the next coin.
//
  while ( 0 < total )
  {
    if ( coin_value[j] <= total )
    {
      total = total - coin_value[j];
      change[change_num] = j;
      change_num = change_num + 1;
    }
    else
    {
      j = j - 1;
      if ( j < 0 )
      {
        break;
      }
    }
  }
  return;
}
//****************************************************************************80

void change_next ( int total, int coin_num, int coin_value[], int &change_num, 
  int change[], bool &done  )

//****************************************************************************80
//
//  Purpose:
//
//    CHANGE_NEXT computes the next set of change for a given sum.
//
//  Example:
//
//    Total = 17
//    COIN_NUM = 3
//    COIN_VALUE = { 1, 5, 10 }
//
//
//        #  CHANGE              COIN_VALUE(CHANGE)
//
//    1   4  3 2 1 1             10 5 1 1
//    2   8  3 1 1 1 1 1 1 1     10 1 1 1 1 1 1 1
//    3   5  2 2 2 1 1            5 5 5 1 1
//    4   9  2 2 1 1 1 1 1 1 1    5 5 1 1 1 1 1 1 1
//    5  13  2 1 1 1 1 1 1 1 1 1  5 1 1 1 1 1 1 1 1 1
//           1 1 1                1 1 1
//    6  17  1 1 1 1 1 1 1 1 1 1  1 1 1 1 1 1 1 1 1 1 1
//           1 1 1 1 1 1 1        1 1 1 1 1 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TOTAL, the total for which change is to be made.
//
//    Input, int COIN_NUM, the number of types of coins.
//
//    Input, int COIN_VALUE[COIN_NUM], the value of each coin.
//    The values must be in ascending order.
//
//    Input/output, int &CHANGE_NUM, the number of coins given in change
//    for this form of the change.
//
//    Input/output, int CHANGE[CHANGE_NUM], the indices of the coins.
//    The user must dimension this array to have dimension TOTAL!
//
//    Input/output, bool &DONE.  The user sets DONE = .TRUE. on
//    first call to tell the routine this is the beginning of a computation.
//    The program resets DONE to .FALSE. and it stays that way until
//    the last possible change combination is made, at which point the
//    program sets DONE to TRUE again.
//
{
  int change_num2;
  int coin_num2;
  int i;
  int last;
  int total2;

  if ( done )
  {
//
//  Make sure the coin values are sorted into ascending order.
//
    if ( !i4vec_ascends ( coin_num, coin_value ) )
    {
      cerr << "\n";
      cerr << "CHANGE_NEXT - Fatal error!\n";
      cerr << "  COIN_VALUE array is not in ascending order.\n";
      exit ( 1 );
    }
//
//  Start with the greedy change.
//
    change_greedy ( total, coin_num, coin_value, change_num, change );
//
//  In a few cases, like change for 4 cents, we're done after the first call.
//
    if ( change_num == total )
    {
      done = true;
    }
    else
    {
      done = false;
    }
    return;
  }
//
//  Find the last location in the input change which is NOT a penny.
//
  last = -1;

  for ( i = change_num-1; 0 <= i; i-- )
  {
    if ( change[i] != 0 )
    {
      last = i;
      break;
    }
  }
//
//  If that location is still -1, an error was made.
//
  if ( last == -1 )
  {
    done = true;
    return;
  }
//
//  Sum the entries from that point to the end.
//
  total2 = 0;
  for ( i = last; i <= change_num-1; i++ )
  {
    total2 = total2 + coin_value [ change[i] ];
  }
//
//  Make greedy change for the partial sum using coins smaller than that one.
//
  coin_num2 = change[last];

  change_greedy ( total2, coin_num2, coin_value, change_num2,
    change+last );

  change_num = last + change_num2;

  return;
}
//****************************************************************************80

bool chinese_check ( int n, int m[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHINESE_CHECK checks the Chinese remainder moduluses.
//
//  Discussion:
//
//    For a Chinese remainder representation, the moduluses M(I) must
//    be positive and pairwise prime.  Also, in case this is not obvious,
//    no more than one of the moduluses may be 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of moduluses.
//
//    Input, int M[N], the moduluses.  These should be positive
//    and pairwise prime.
//
//    Output, bool CHINESE_CHECK, is TRUE if an error was detected.
//
{
  int i;
  int j;
//
//  Do not allow nonpositive entries.
//
  for ( i = 0; i < n; i++ )
  {
    if ( m[i] <= 0 )
    {
      return true;
    }
  }
//
//  Allow one entry to be 1, but not two entries.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = i+1; j < n; j++ )
    {
      if ( m[i] == 1 && m[j] == 1 )
      {
        return true;
      }
    }
  }
//
//  Now check pairwise primeness.
//
  if ( !i4vec_pairwise_prime ( n, m ) )
  {
    return true;
  }

  return false;
}
//****************************************************************************80

int chinese_to_i4 ( int n, int m[], int r[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHINESE_TO_I4 converts a set of Chinese remainders to an equivalent integer.
//
//  Discussion:
//
//    Given a set of N pairwise prime, positive moduluses M(I), and
//    a corresponding set of remainders R(I), this routine finds an
//    integer J such that, for all I,
//
//      J = R(I) mod M(I)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of moduluses.
//
//    Input, int M[N], the moduluses.  These should be positive
//    and pairwise prime.
//
//    Input, int R[N], the Chinese remainder representation of the integer.
//
//    Output, int CHINESE_TO_I4, the corresponding integer.
//
{
  int a;
  int *b;
  int big_m;
  int c;
  bool error;
  int i;
  int j;

  error = chinese_check ( n, m );

  if ( error )
  {
    cerr << "\n";
    cerr << "CHINESE_TO_I4 - Fatal error!\n";
    cerr << "  The moduluses are not legal.\n";
    exit ( 1 );
  }

  b = new int[n];
//
//  Set BIG_M.
//
  big_m = 1;
  for ( i = 0; i < n; i++ )
  {
    big_m = big_m * m[i];
  }
//
//  Solve BIG_M / M(I) * B(I) = 1, mod M(I)
//
  for ( i = 0; i < n; i++ )
  {
    a = big_m / m[i];
    c = 1;
    b[i] = congruence ( a, m[i], c, error );
  }
//
//  Set J = sum ( 1 <= I <= N ) ( R(I) * B(I) * BIG_M / M(I) ) mod M
//
  j = 0;
  for ( i = 0; i < n; i++ )
  {
    j = ( j + r[i] * b[i] * ( big_m / m[i] ) ) % big_m;
  }

  delete [] b;

  return j;
}
//****************************************************************************80

void comb_next ( int n, int k, int a[], bool &done )

//****************************************************************************80
//
//  Purpose:
//
//    COMB_NEXT computes combinations of K things out of N.
//
//  Discussion:
//
//    The combinations are computed one at a time, in lexicographical order.
//
//    10 April 1009: Thanks to "edA-qa mort-ora-y" for supplying a 
//    correction to this code!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Mifsud,
//    Combination in Lexicographic Order,
//    ACM algorithm 154,
//    Communications of the ACM,
//    March 1963.
//
//  Parameters:
//
//    Input, int N, the total number of things.
//
//    Input, int K, the number of things in each combination.
//
//    Input/output, int A[K], contains the list of elements in
//    the current combination.
//
//    Input/output, bool &DONE.  Set DONE to TRUE before the first call,
//    and then use the output value from the previous call on subsequent
//    calls.  The output value will be FALSE as long as there are more
//    combinations to compute, and TRUE when the list is exhausted.
//
{
  int i;
  int j;

  if ( done )
  {
    if ( k <= 0 )
    {
      return;
    }

    i4vec_indicator ( k, a );

    done = false; 
  }
  else
  {
    if ( a[k-1] < n )
    {
      a[k-1] = a[k-1] + 1;
      return;
    }

    for ( i = k; 2 <= i; i-- )
    {
      if ( a[i-2] < n-k+i-1 )
      {
        a[i-2] = a[i-2] + 1;

        for ( j = i; j <= k; j++ )
        {
          a[j-1] = a[i-2] + j - ( i-1 );
        }
        return;
      }
    }
    done = true;
  }

  return;
}
//****************************************************************************80

void comb_row ( bool next, int n, int row[] )

//****************************************************************************80
//
//  Purpose:
//
//    COMB_ROW computes row N of Pascal's triangle.
//
//  Discussion:
//
//    Row N contains the N+1 combinatorial coefficients
//
//      C(N,0), C(N,1), C(N,2), ... C(N,N)
//
//    The sum of the elements of row N is equal to 2**N.
//
//    The formula is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  First terms:
//
//     N K:0  1   2   3   4   5   6   7  8  9 10
//
//     0   1
//     1   1  1
//     2   1  2   1
//     3   1  3   3   1
//     4   1  4   6   4   1
//     5   1  5  10  10   5   1
//     6   1  6  15  20  15   6   1
//     7   1  7  21  35  35  21   7   1
//     8   1  8  28  56  70  56  28   8  1
//     9   1  9  36  84 126 126  84  36  9  1
//    10   1 10  45 120 210 252 210 120 45 10  1
//
//  Recursion:
//
//    C(N,K) = C(N-1,K-1)+C(N-1,K)
//
//  Special values:
//
//    C(N,0) = C(N,N) = 1
//    C(N,1) = C(N,N-1) = N
//    C(N,N-2) = sum ( 1 <= I <= N ) N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, bool NEXT, indicates whether this is a call for
//    the 'next' row of the triangle.
//    NEXT = FALSE means this is a startup call.  Row N is desired, but
//    presumably this is a first call, or row N-1 was not computed
//    on the previous call.
//    NEXT = TRUE means this is not the first call, and row N-1 was computed
//    on the previous call.  In this case, much work can be saved
//    by using the information from the previous values of ROW
//    to build the next values.
//
//    Input, int N, the row of the triangle desired.  The triangle
//    begins with row 0.
//
//    Output, int ROW[N+1], the row of coefficients.
//    ROW(I) = C(N,I-1).
//
{
  int i;
  int j;

  if ( n < 0 )
  {
    return;
  }

  if ( next )
  {
    for ( i = n-1; 1 <= i; i-- )
    {
      row[i] = row[i] + row[i-1];
    }
    row[n] = 1;
  }
  else
  {
    row[0] = 1;
    for ( i = 1; i <= n; i++ )
    {
      row[i] = 0;
    }

    for ( j = 1; j <= n; j++ )
    {
      for ( i = j; 1 <= i; i-- )
      {
        row[i] = row[i] + row[i-1];
      }
    }
  }

  return;
}
//****************************************************************************80

void comb_unrank ( int m, int n, int rank, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    COMB_UNRANK returns the RANK-th combination of N things out of M.
//
//  Discussion:
//
//    The combinations are ordered lexically.
//
//    Lexical order can be illustrated for the general case of N and M as
//    follows:
//
//    1:       1,     2,     3,     ..., N-2, N-1, N
//    2:       1,     2,     3,     ..., N-2, N-1, N+1
//    3:       1,     2,     3,     ..., N-2, N-1, N+2
//    ...
//    M-N+1:   1,     2,     3,     ..., N-2, N-1, M
//    M-N+2:   1,     2,     3,     ..., N-2, N,   N+1
//    M-N+3:   1,     2,     3,     ..., N-2, N,   N+2
//    ...
//    LAST-2:  M-N,   M-N+1, M-N+3, ..., M-2, M-1, M
//    LAST-1:  M-N,   M-N+2, M-N+3, ..., M-2, M-1, M
//    LAST:    M-N+1, M-N+2, M-N+3, ..., M-2, M-1, M
//
//    There are a total of M!/(N!*(M-N)!) combinations of M
//    things taken N at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    B P Buckles, M Lybanon,
//    Algorithm 515,
//    Generation of a Vector from the Lexicographical Index,
//    ACM Transactions on Mathematical Software,
//    Volume 3, Number 2, pages 180-182, June 1977.
//
//  Parameters:
//
//    Input, int M, the size of the set.
//
//    Input, int N, the number of things in the combination.
//    N must be greater than 0, and no greater than M.
//
//    Input, int RANK, the lexicographical index of combination
//    sought.  RANK must be at least 1, and no greater than M!/(N!*(M-N)!).
//
//    Output, int A[N], array containing the combination set.
//
{
  int i;
  int j;
  int k;
//
//  Initialize lower bound index at zero.
//
  k = 0;
//
//  Loop to select elements in ascending order.
//
  for ( i = 1; i <= n-1; i++ )
  {
//
//  Set lower bound element number for next element value.
//
    a[i-1] = 0;

    if ( 1 < i )
    {
      a[i-1] = a[i-2];
    }
//
//  Check each element value.
//
    for ( ; ; )
    {
      a[i-1] = a[i-1] + 1;
      j = i4_choose ( m-a[i-1], n-i );
      k = k + j;

      if ( rank <= k )
      {
        break;
      }
    }
    k = k - j;
  }

  a[n-1] = a[n-2] + rank - k;

  return;
}
//****************************************************************************80

void comp_next ( int n, int k, int a[], bool &more, int &h, int &t )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT computes the compositions of the integer N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is
//    3+2+1, another would be 6+0+0.
//
//    On the first call to this routine, set MORE = FALSE.  The routine
//    will compute the first element in the sequence of compositions, and
//    return it, as well as setting MORE = TRUE.  If more compositions
//    are desired, call again, and again.  Each time, the routine will
//    return with a new composition.
//
//    However, when the LAST composition in the sequence is computed 
//    and returned, the routine will reset MORE to FALSE, signaling that
//    the end of the sequence has been reached.
//
//    This routine originally used a SAVE statement to maintain the
//    variables H and T.  I have decided that it is safer
//    to pass these variables as arguments, even though the user should
//    never alter them.  This allows this routine to safely shuffle
//    between several ongoing calculations.
//
//
//    There are 28 compositions of 6 into three parts.  This routine will
//    produce those compositions in the following order:
//
//     I         A
//     -     ---------
//     1     6   0   0
//     2     5   1   0
//     3     4   2   0
//     4     3   3   0
//     5     2   4   0
//     6     1   5   0
//     7     0   6   0
//     8     5   0   1
//     9     4   1   1
//    10     3   2   1
//    11     2   3   1
//    12     1   4   1
//    13     0   5   1
//    14     4   0   2
//    15     3   1   2
//    16     2   2   2
//    17     1   3   2
//    18     0   4   2
//    19     3   0   3
//    20     2   1   3
//    21     1   2   3
//    22     0   3   3
//    23     2   0   4
//    24     1   1   4
//    25     0   2   4
//    26     1   0   5
//    27     0   1   5
//    28     0   0   6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2008
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool &MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
//    Input/output, int &H, &T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;

  if ( !( more ) )
  {
    t = n;
    h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < t )
    {
      h = 0;
    }
    h = h + 1;
    t = a[h-1];
    a[h-1] = 0;
    a[0] = t - 1;
    a[h] = a[h] + 1;
  }

  more = ( a[k-1] != n );

  return;
}
//****************************************************************************80

void comp_random ( int n, int k, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_RANDOM selects a random composition of the integer N into K parts.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer to be decomposed.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int A[K], the parts of the composition.
//
{
  int i;
  int l;
  int m;

  ksub_random ( n+k-1, k-1, seed, a );

  a[k-1] = n + k;
  l = 0;

  for ( i = 0; i < k; i++ )
  {
    m = a[i];
    a[i] = a[i] - l - 1;
    l = m;
  }

  return;
}
//****************************************************************************80

void compnz_next ( int n, int k, int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    COMPNZ_NEXT computes the compositions of the integer N into K nonzero parts.
//
//  Discussion:
//
//    A composition of the integer N into K nonzero parts is an ordered sequence
//    of K positive integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is 3+2+1, another would
//    be 4+1+1 but 5+1+0 is not allowed since it includes a zero part.
//
//    On the first call to this routine, set MORE = FALSE.  The routine
//    will compute the first element in the sequence of compositions, and
//    return it, as well as setting MORE = TRUE.  If more compositions
//    are desired, call again, and again.  Each time, the routine will
//    return with a new composition.
//
//    However, when the LAST composition in the sequence is computed
//    and returned, the routine will reset MORE to FALSE, signaling that
//    the end of the sequence has been reached.
//
//  Example:
//
//    The 10 compositions of 6 into three nonzero parts are:
//
//      4 1 1,  3 2 1,  3 1 2,  2 3 1,  2 2 2,  2 1 3,
//      1 4 1,  1 3 2,  1 2 3,  1 1 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//    K must be no greater than N.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool &MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
{
  int i;
  static int h = 0;
  static int t = 0;
//
//  We use the trick of computing ordinary compositions of (N-K)
//  into K parts, and adding 1 to each part.
//
  if ( n < k )
  {
    more = false;
    for ( i = 0; i < k; i++ )
    {
      a[i] = -1;
    }
    return;
  }
//
//  The first computation.
//
  if ( !( more ) )
  {
    t = n - k;
    h = 0;
    a[0] = n - k;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    for ( i = 0; i < k; i++ )
    {
      a[i] = a[i] - 1;
    }
    if ( 1 < t )
    {
      h = 0;
    }

    h = h + 1;
    t = a[h-1];
    a[h-1] = 0;
    a[0] = t - 1;
    a[h] = a[h] + 1;
  }

  more = ( a[k-1] != ( n - k ) );

  for ( i = 0; i < k; i++ )
  {
    a[i] = a[i] + 1;
  }

  return;
}
//****************************************************************************80

void compnz_random ( int n, int k, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    COMPNZ_RANDOM selects a random composition of the integer N into K nonzero parts.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer to be decomposed.
//
//    Input, int K, the number of parts in the composition.
//    K must be no greater than N.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int A[K], the parts of the composition.
//
{
  int i;
  int l;
  int m;

  if ( n < k )
  {
    for ( i = 0; i < k; i++ )
    {
      a[i] = -1;
    }
    return;
  }

  ksub_random ( n-1, k-1, seed, a );

  a[k-1] = n;
  l = 0;

  for ( i = 0; i < k; i++ )
  {
    m = a[i];
    a[i] = a[i] - l - 1;
    l = m;
  }

  for ( i = 0; i < k; i++ )
  {
    a[i] = a[i] + 1;
  }

  return;
}
//****************************************************************************80

int congruence ( int a, int b, int c, bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    CONGRUENCE solves a congruence of the form ( A * X = C ) mod B.
//
//  Discussion:
//
//    A, B and C are given integers.  The equation is solvable if and only
//    if the greatest common divisor of A and B also divides C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric Weisstein, editor,
//    CRC Concise Encylopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Parameters:
//
//    Input, int A, B, C, the coefficients of the Diophantine equation.
//
//    Output, bool &ERROR, error flag, is TRUE if an error occurred..
//
//    Output, int CONGRUENCE, the solution of the Diophantine equation.
//    X will be between 0 and B-1.
//
{
# define N_MAX 100

  int a_copy;
  int a_mag;
  int a_sign;
  int b_copy;
  int b_mag;
  int b_sign;
  int c_copy;
  int g;
  int k;
  int n;
  float norm_new;
  float norm_old;
  int q[N_MAX];
  bool swap;
  int temp;
  int x;
  int xnew;
  int y;
  int ynew;
  int z;
//
//  Defaults for output parameters.
//
  error = false;
  x = 0;
  y = 0;
//
//  Special cases.
//
  if ( a == 0 && b == 0 && c == 0 )
  {
    x = 0;
    return x;
  }
  else if ( a == 0 && b == 0 && c != 0 )
  {
    error = true;
    x = 0;
    return x;
  }
  else if ( a == 0 && b != 0 && c == 0 )
  {
    x = 0;
    return x;
  }
  else if ( a == 0 && b != 0 && c != 0 )
  {
    x = 0;
    if ( ( c % b ) != 0 )
    {
      error = true;
    }
    return x;
  }
  else if ( a != 0 && b == 0 && c == 0 )
  {
    x = 0;
    return x;
  }
  else if ( a != 0 && b == 0 && c != 0 )
  {
    x = c / a;
    if ( ( c % a ) != 0 )
    {
      error = true;
    }
    return x;
  }
  else if ( a != 0 && b != 0 && c == 0 )
  {
//  g = i4_gcd ( a, b );
//  x = b / g;
    x = 0;
    return x;
  }
//
//  Now handle the "general" case: A, B and C are nonzero.
//
//  Step 1: Compute the GCD of A and B, which must also divide C.
//
  g = i4_gcd ( a, b );

  if ( ( c % g ) != 0 )
  {
    error = true;
    return x;
  }

  a_copy = a / g;
  b_copy = b / g;
  c_copy = c / g;
//
//  Step 2: Split A and B into sign and magnitude.
//
  a_mag = abs ( a_copy );
  a_sign = i4_sign ( a_copy );
  b_mag = abs ( b_copy );
  b_sign = i4_sign ( b_copy );
//
//  Another special case, A_MAG = 1 or B_MAG = 1.
//
  if ( a_mag == 1 )
  {
    x = a_sign * c_copy;
    return x;
  }
  else if ( b_mag == 1 )
  {
    x = 0;
    return x;
  }
//
//  Step 3: Produce the Euclidean remainder sequence.
//
  if ( b_mag <= a_mag )
  {
    swap = false;
    q[0] = a_mag;
    q[1] = b_mag;
  }
  else
  {
    swap = true;
    q[0] = b_mag;
    q[1] = a_mag;
  }

  n = 3;

  for ( ; ; )
  {
    q[n-1] = ( q[n-3] % q[n-2] );

    if ( q[n-1] == 1 )
    {
      break;
    }

    n = n + 1;

    if ( N_MAX < n )
    {
      error = true;
      cerr << "\n";
      cerr << "CONGRUENCE - Fatal error!\n";
      cerr << "  Exceeded number of iterations.\n";
      exit ( 1 );
    }
  }
//
//  Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
//
  y = 0;
  for ( k = n; 2 <= k; k-- )
  {
    x = y;
    y = ( 1 - x * q[k-2] ) / q[k-1];
  }
//
//  Step 5: Undo the swapping.
//
  if ( swap )
  {
    z = x;
    x = y;
    y = z;
  }
//
//  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
//
  x = x * a_sign;
//
//  Step 7: Multiply by C, so that X * A + Y * B = C.
//
  x = x * c_copy;
//
//  Step 8: Now force 0 <= X < B.
//
  x = x % b;
//
//  Step 9: Force positivity.
//
  if ( x < 0 )
  {
    x = x + b;
  }

  return x;
# undef N_MAX
}
//****************************************************************************80

void count_pose_random ( int &seed, int blocks[], int &goal )

//****************************************************************************80
//
//  Purpose:
//
//    COUNT_POSE_RANDOM poses a problem for the game "The Count is Good"
//
//  Discussion:
//
//    The French television show "The Count is Good" has a game that goes
//    as follows:
//
//      A number is chosen at random between 100 and 999.  This is the GOAL.
//
//      Six numbers are randomly chosen from the set 1, 2, 3, 4, 5, 6, 7, 8,
//      9, 10, 25, 50, 75, 100.  These numbers are the BLOCKS.
//
//      The player must construct a formula, using some or all of the blocks,
//      (but not more than once), and the operations of addition, subtraction,
//      multiplication and division.  Parentheses should be used to remove
//      all ambiguity.  However, it is forbidden to use subtraction in a
//      way that produces a negative result, and all division must come out
//      exactly, with no remainder.
//
//    This routine poses a sample problem from the show.  The point is,
//    to determine how to write a program that can solve such a problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Raymond Seroul,
//    Programming for Mathematicians,
//    Springer Verlag, 2000, page 355-357.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int BLOCKS[6], the six numbers available for the formula.
//
//    Output, int &GOAL, the goal number.
//
{
  int i;
  int ind[6];
  static int stuff[14] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 75, 100 };

  goal = i4_uniform ( 100, 999, seed );

  ksub_random ( 14, 6, seed, ind );

  for ( i = 0; i < 6; i++ )
  {
    blocks[i] = stuff[ind[i]-1];
  }

  return;
}
//****************************************************************************80

void debruijn ( int m, int n, int string[] )

//****************************************************************************80
//
//  Purpose:
//
//    DEBRUIJN constructs a de Bruijn sequence.
//
//  Discussion:
//
//    Suppose we have an alphabet of M letters, and we are interested in
//    all possible strings of length N.  If M = 2 and N = 3, then we are
//    interested in the M**N strings:
//
//      000
//      001
//      010
//      011
//      100
//      101
//      110
//      111
//
//    Now, instead of making a list like this, we prefer, if possible, to
//    write a string of letters, such that every consecutive sequence of
//    N letters is one of the strings, and every string occurs once, if
//    we allow wraparound.
//
//    For the above example, a suitable sequence would be the 8 characters:
//
//      00011101(00...
//
//    where we have suggested the wraparound feature by repeating the first
//    two characters at the end.
//
//    Such a sequence is called a de Bruijn sequence.  It can easily be
//    constructed by considering a directed graph, whose nodes are all
//    M**(N-1) strings of length N-1.  A node I has a directed edge to
//    node J (labeled with character K) if the string at node J can
//    be constructed by beheading the string at node I and adding character K.
//
//    In this setting, a de Bruijn sequence is simply an Eulerian circuit
//    of the directed graph, with the edge labels being the entries of the
//    sequence.  In general, there are many distinct de Bruijn sequences
//    for the same parameter M and N.  This program will only find one
//    of them.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of letters in the alphabet.
//
//    Input, int N, the number of letters in a codeword.
//
//    Output, int STRING[M**N], a deBruijn string.
//
{
  int i;
  int iedge;
  int *inode;
  int *ivec;
  int j;
  int *jnode;
  int *jvec;
  int k;
  int *knode;
  int nedge;
  int nnode;
  bool success;
  int *trail;
//
//  Construct the adjacency information.
//
  nnode = i4_power ( m, n-1 );
  nedge = i4_power ( m, n );

  inode = new int[nedge];
  ivec = new int[n-1];
  jnode = new int[nedge];
  jvec = new int[n-1];
  knode = new int[nedge];

  iedge = 0;

  for ( i = 1; i <= nnode; i++ )
  {
    index_unrank0 ( n-1, m, i, ivec );

    for ( k = 1; k <= m; k++ )
    {
//
//  Shift N-2 entries of IVEC down.
//
      for ( j = 0; j < n-2; j++ )
      {
        jvec[j] = ivec[j+1];
      }
      jvec[n-2] = k;

      j = index_rank0 ( n-1, m, jvec );

      inode[iedge] = i;
      jnode[iedge] = j;
      knode[iedge] = k;
      iedge = iedge + 1;
    }
  }

  delete [] ivec;
  delete [] jvec;
//
//  Determine a circuit.
//
  trail = new int[nedge];

  digraph_arc_euler ( nnode, nedge, inode, jnode, success, trail );
//
//  The string is constructed from the labels of the edges in the circuit.
//
  for ( i = 0; i < nedge; i++ )
  {
    string[i] = knode[trail[i]-1];
  }

  delete [] inode;
  delete [] jnode;
  delete [] knode;
  delete [] trail;

  return;
}
//****************************************************************************80

void dec_add ( int mantissa1, int exponent1, int mantissa2, int exponent2, 
  int dec_digit, int &mantissa, int &exponent )

//****************************************************************************80
//
//  Purpose:
//
//    DEC_ADD adds two decimal quantities.
//
//  Discussion:
//
//    A decimal value is represented as MANTISSA * 10^EXPONENT.
//
//    The routine computes
//
//      MANTISSA * 10^EXPONENT = MANTISSA1 * 10^EXPONENT1 + MANTISSA2 * 10^EXPONENT2
//
//    using DEC_DIGIT arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MANTISSA1, EXPONENT1, the first number to be added.
//
//    Input, int MANTISSA2, EXPONENT2, the second number to be added.
//
//    Input, int DEC_DIGIT, the number of decimal digits.
//
//    Output, int &MANTISSA, &EXPONENT, the sum.
//
{
  if ( mantissa1 == 0 )
  {
    mantissa = mantissa2;
    exponent = exponent2;
    dec_round ( mantissa, exponent, dec_digit, mantissa, exponent );
    return;
  }
  else if ( mantissa2 == 0 )
  {
    mantissa = mantissa1;
    exponent = exponent1;
    dec_round ( mantissa, exponent, dec_digit, mantissa, exponent );
    return;
  }
//
//  Line up the exponents.
//
  if ( exponent1 < exponent2 )
  {
    mantissa2 = mantissa2 * ( int ) pow ( ( double ) 10, ( exponent2 - exponent1 ) );
    exponent2 = exponent1;
    mantissa = mantissa1 + mantissa2;
    exponent = exponent1;
  }
  else if ( exponent1 == exponent2 ) 
  {
    mantissa = mantissa1 + mantissa2;
    exponent = exponent1;
  }
  else if ( exponent2 < exponent1 )
  {
    mantissa1 = mantissa1 * ( int ) pow ( ( double ) 10, ( exponent1 - exponent2 ) );
    exponent1 = exponent2;
    mantissa = mantissa1 + mantissa2;
    exponent = exponent2;
  }
//
//  Clean up the result.
//
  dec_round ( mantissa, exponent, dec_digit, mantissa, exponent );

  return;
}
//****************************************************************************80

void dec_div ( int mantissa1, int exponent1, int mantissa2, int exponent2,
  int dec_digit, int &mantissa, int &exponent, bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    DEC_DIV divides two decimal values.
//
//  Discussion:
//
//    A decimal value is represented as MANTISSA * 10^EXPONENT.
//
//    The routine computes
//
//      MANTISSA * 10**EXPONENT 
//      = (MANTISSA1 * 10**EXPONENT1) / (MANTISSA2 * 10**EXPONENT2)
//      = (MANTISSA1/MANTISSA2) * 10**(EXPONENT1-EXPONENT2)
//
//    while avoiding integer overflow.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MANTISSA1, EXPONENT1, the numerator.
//
//    Input, int MANTISSA2, EXPONENT2, the denominator.
//
//    Input, int DEC_DIGIT, the number of decimal digits.
//
//    Output, int &MANTISSA, &EXPONENT, the result.
//
//    Output, bool &ERROR is true if an error occurred.
//
{
  double dval;
  int exponent3;
  int mantissa3;

  error = false;
//
//  First special case, top fraction is 0.
//
  if ( mantissa1 == 0 )
  {
    mantissa = 0;
    exponent = 0;
    return;
  }
//
//  First error, bottom of fraction is 0.
//
  if ( mantissa2 == 0 )
  {
    error = true;
    mantissa = 0;
    exponent = 0;
    return;
  }
//
//  Second special case, result is 1.
//
  if ( mantissa1 == mantissa2 && exponent1 == exponent2 )
  {
    mantissa = 1;
    exponent = 0;
    return;
  }
//
//  Third special case, result is power of 10.
//
  if ( mantissa1 == mantissa2 )
  {
    mantissa = 1;
    exponent = exponent1 - exponent2;
    return;
  }
//
//  Fourth special case: MANTISSA1/MANTISSA2 is exact.
//
  if ( ( mantissa1 / mantissa2 ) * mantissa2 == mantissa1 )
  {
    mantissa = mantissa1 / mantissa2;
    exponent = exponent1 - exponent2;
    return;
  }
//
//  General case.
//
  dval = ( double ) mantissa1 / ( double ) mantissa2;

  r8_to_dec ( dval, dec_digit, mantissa3, exponent3 );

  mantissa = mantissa3;
  exponent = exponent3 + exponent1 - exponent2;

  return;
}
//****************************************************************************80

void dec_mul ( int mantissa1, int exponent1, int mantissa2, int exponent2, 
  int dec_digit, int &mantissa, int &exponent )

//****************************************************************************80
//
//  Purpose:
//
//    DEC_MUL multiplies two decimals.
//
//  Discussion:
//
//    A decimal value is represented as MANTISSA * 10**EXPONENT.
//
//    The routine computes
//
//      MANTISSA * 10**EXPONENT = (MANTISSA1 * 10**EXPONENT1) * (MANTISSA2 * 10**EXPONENT2)
//                      = (MANTISSA1*MANTISSA2) * 10**(EXPONENT1+EXPONENT2)
//
//    while avoiding integer overflow.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MANTISSA1, EXPONENT1, the first multiplier.
//
//    Input, int MANTISSA2, EXPONENT2, the second multiplier.
//
//    Input, int DEC_DIGIT, the number of decimal digits.
//
//    Output, int &MANTISSA, &EXPONENT, the product.
//
{
  double dval;
  int exponent3;
  int mantissa3;
  double temp;
//
//  The result is zero if either MANTISSA1 or MANTISSA2 is zero.
//
  if ( mantissa1 == 0 || mantissa2 == 0 )
  {
    mantissa = 0;
    exponent = 0;
    return;
  }
//
//  The result is simple if either MANTISSA1 or MANTISSA2 is one.
//
  if ( abs ( mantissa1 ) == 1 || abs ( mantissa2 ) == 1 )
  {
    mantissa = mantissa1 * mantissa2;
    exponent = exponent1 + exponent2;
    return;
  }

  temp = log ( ( double ) abs ( mantissa1 ) ) 
       + log ( ( double ) abs ( mantissa2 ) );

  if ( temp < log ( ( double ) i4_huge ( ) ) )
  {
    mantissa = mantissa1 * mantissa2;
    exponent = exponent1 + exponent2;
  }
  else
  {
    dval = ( double ) mantissa1 * ( double ) mantissa2;

    r8_to_dec ( dval, dec_digit, mantissa3, exponent3 );

    mantissa = mantissa3;
    exponent = exponent3 + ( exponent1 + exponent2 );
  }

  dec_round ( mantissa, exponent, dec_digit, mantissa, exponent );

  return;
}
//****************************************************************************80

void dec_round ( int mantissa1, int exponent1, int dec_digit, 
  int &mantissa2, int &exponent2 )

//****************************************************************************80
//
//  Purpose:
//
//    DEC_ROUND rounds a decimal fraction to a given number of digits.
//
//  Discussion:
//
//    A decimal value is represented as MANTISSA * 10^EXPONENT.
//
//    The routine takes an arbitrary decimal value and makes sure that MANTISSA
//    has no more than DEC_DIGIT digits.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MANTISSA1, *EXPONENT1, the coefficient and exponent
//    of a decimal fraction to be rounded.
//
//    Input, int DEC_DIGIT, the number of decimal digits.
//
//    Output, int &MANTISSA2, &EXPONENT2, the rounded coefficient and exponent
//    of a decimal fraction.  MANTISSA2 has no more than
//    DEC_DIGIT decimal digits.
//
{
  int i;
  int limit;
  int sgn;

  mantissa2 = mantissa1;
  exponent2 = exponent1;
//
//  Watch out for the special case of 0.
//
  if ( mantissa2 == 0 )
  {
    exponent2 = 0;
    return;
  }
//
//  Record the sign of MANTISSA.
//
  sgn = 1;
  if ( mantissa2 < 0 ) 
  {
    mantissa2 = -( mantissa2 );
    sgn = -sgn;
  }
//
//  If MANTISSA is too big, knock it down.
//
  limit = 1;
  for ( i = 1; i <= dec_digit; i++ )
  {
    limit = limit * 10;
  }

  while ( limit <= abs ( mantissa2 ) )
  {
    mantissa2 = ( mantissa2 + 5 ) / 10;
    exponent2 = exponent2 + 1;
  }
//
//  Absorb trailing 0's into the exponent.
//
  while ( ( mantissa2 / 10 ) * 10 == mantissa2 )
  {
    mantissa2 = mantissa2 / 10;
    exponent2 = exponent2 + 1;
  }

  mantissa2 = sgn * ( mantissa2 );

  return;
}
//****************************************************************************80

double dec_to_r8 ( int mantissa, int exponent )

//****************************************************************************80
//
//  Purpose:
//
//    DEC_TO_R8 converts a decimal value to an R8.
//
//  Discussion:
//
//    A decimal value is represented as MANTISSA * 10**EXPONENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MANTISSA, EXPONENT, the coefficient and exponent
//    of the decimal value.
//
//    Output, double DEC_TO_R8, the real value of the decimal.
//
{
  double value;

  value = mantissa * pow ( 10.0, exponent );

  return value;
}
//****************************************************************************80

void dec_to_rat ( int mantissa, int exponent, int &rat_top, int &rat_bot )

//****************************************************************************80
//
//  Purpose:
//
//    DEC_TO_RAT converts a decimal to a rational representation.
//
//  Discussion:
//
//    A decimal value is represented as MANTISSA * 10**EXPONENT.
//
//    A rational value is represented by RAT_TOP / RAT_BOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MANTISSA, EXPONENT, the decimal number.
//
//    Output, int &RAT_TOP, &RAT_BOT, the rational value.
//
{
  int gcd;
  int i;

  if ( exponent == 0 )
  {
    rat_top = mantissa;
    rat_bot = 1;
  }
  else if ( 0 < exponent )
  {
    rat_top = mantissa;
    for ( i = 1; i <= exponent; i++ )
    {
      rat_top = rat_top * 10;
    }
    rat_bot = 1;
  }
  else
  {
    rat_top = mantissa;
    rat_bot = 1;
    for ( i = 1; i <= -exponent; i++ )
    {
      rat_bot = rat_bot * 10;
    }

    gcd = i4_gcd ( rat_top, rat_bot );
    rat_top = rat_top / gcd;
    rat_bot = rat_bot / gcd;
  }

  return;
}
//****************************************************************************80

char *dec_to_s ( int mantissa, int exponent )

//****************************************************************************80
//
//  Purpose:
//
//    DEC_TO_S converts a decimal number to a string.
//
//  Discussion:
//
//    This is a draft version that is NOT WORKING YET.
//
//  Example:
//
//    Mantissa  Exponent Representation:
//
//         523      -1              5.23
//         134       2          13400
//           0      10              0
//
//      123456       3      123456000
//      123456       2       12345600
//      123456       1        1234560
//      123456       0         123456
//      123456      -1          12345.6
//      123456      -2           1234.56
//      123456      -3            123.456
//      123456      -4             12.3456
//      123456      -5             1.23456
//      123456      -6             0.123456
//      123456      -7             0.0123456
//      123456      -8             0.00123456
//      123456      -9             0.000123456
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MANTISSA, EXPONENT, integers which represent the decimal.
//
//    Output, char *S, the representation of the value.
//
{
  int digit;
  int i;
  int last;
  int mantissa_exponent;
  int mantissa_exponent_copy;
  int mantissa_10;
  int pos;
  char *s;
  int s_length;
  int sign;

  s_length = dec_width ( mantissa, exponent ) + 1;

  s = new char[s_length];

  for ( i = 0; i < s_length - 1; i++ )
  {
    s[i] = '0';
  }
  s[s_length-1] = '\0';

  if ( mantissa == 0 )
  {
    return s;
  }

  pos = 0;

  if ( mantissa < 0 )
  {
    s[pos] = '-';
    pos = pos + 1;
    mantissa = -mantissa;
  }

  mantissa_exponent = i4_log_10 ( mantissa ) + 1;
  mantissa_10 = ( int ) pow ( ( double ) 10, mantissa_exponent - 1 );
//
//  Are the next characters "0."?
//
  if ( mantissa_exponent + exponent <= 0 ) 
  {
    s[pos] = '0';
    pos = pos + 1;
    s[pos] = '.';
    pos = pos + 1;

    for ( i = mantissa_exponent + exponent; i < 0; i++ )
    {
      s[pos] = '0';
      pos = pos + 1;
    }
  }
//
//  Print the digits of the mantissa.
//
  mantissa_exponent_copy = mantissa_exponent;

  for ( i = 0; i < mantissa_exponent; i++ )
  {
    digit = mantissa / mantissa_10;
    mantissa = mantissa % mantissa_10;
    s[pos] = digit_to_ch ( digit );
    pos = pos + 1;
    mantissa_10 = mantissa_10 / 10;
    mantissa_exponent_copy = mantissa_exponent_copy - 1;
    if ( exponent < 0 )
    {
      if ( mantissa_exponent_copy + exponent == 0 )
      {
        s[pos] = '.';
        pos = pos + 1;
      }
    }
  }
//
//  Print any trailing zeros.
//

  if ( 0 < exponent )
  {
    for ( i = exponent; 0 < i; i-- )
    {
      s[pos] = '0';
      pos = pos + 1;
    }
  }

  return s;
}
//****************************************************************************80

int dec_width ( int mantissa, int exponent )

//****************************************************************************80
//
//  Purpose:
//
//    DEC_WIDTH returns the "width" of a decimal number.
//
//  Discussion:
//
//    A decimal value is represented as MANTISSA * 10**EXPONENT.
//
//    The "width" of a decimal number is the number of characters
//    required to print it.
//
//  Example:
//
//    Mantissa  Exponent Width  Representation:
//
//         523      -1       4           5.23
//         134       2       5       13400
//           0      10       1           0
//
//      123456       3       9    123456000
//      123456       2       8     12345600
//      123456       1       7      1234560
//      123456       0       6       123456
//      123456      -1       7        12345.6
//      123456      -2       7         1234.56
//      123456      -3       7          123.456
//      123456      -4       7           12.3456
//      123456      -5       7           1.23456
//      123456      -6       8           0.123456
//      123456      -7       9           0.0123456
//      123456      -8      10           0.00123456
//      123456      -9      11           0.000123456
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MANTISSA, EXPONENT, the decimal number.
//
//    Output, int DEC_WIDTH, the "width" of the decimal number.
//
{
  int mantissa_abs;
  int ten_pow;
  int value;
//
//  Special case of 0.
//
  value = 1;

  if ( mantissa == 0 )
  {
    return value;
  }
//
//  Determine a power of 10 that is strictly bigger than MANTISSA.
//  The exponent of that power of 10 is our first estimate for 
//  the number of places.
//
  ten_pow = 10;
  mantissa_abs = abs ( mantissa );

  while ( ten_pow <= mantissa_abs )
  {
    value = value + 1;
    ten_pow = ten_pow * 10;
  }
//
//  If the exponent is nonnegative, that just adds more places.
//
  if ( 0 <= exponent )
  {
    value = value + exponent;
  }
//
//  If the exponent is a little negative, then we are essentially
//  just inserting a decimal point, and moving it around.
//
  else if ( -value < exponent )
  {
    value = value + 1;
  }
//
//  A very negative value of B means we have a leading 0 and decimal,
//  and B trailing places.
//
  else if ( exponent <= -value )
  {
    value = 2 - exponent;
  }
//
//  Take care of sign of MANTISSA.
//
  if ( mantissa < 0 )
  {
    value = value + 1;
  }

  return value;
}
//****************************************************************************80

void decmat_det ( int n, int atop[], int abot[], int dec_digit, 
  int &dtop, int &dbot )

//****************************************************************************80
//
//  Purpose:
//
//    DECMAT_DET finds the determinant of an N by N matrix of decimal entries.
//
//  Discussion:
//
//    The brute force method is used.  The routine should only be used for
//    small matrices, since this calculation requires the summation of N!
//    products of N numbers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of A.
//
//    Input, int ATOP[N*N], ABOT[N*N], the decimal
//    representation of the matrix.
//
//    Input, int DEC_DIGIT, the number of decimal digits.
//
//    Output, int &DTOP, &DBOT, the decimal determinant of the matrix.
//
{
  bool even;
  int i;
  int *iarray;
  int ibot;
  int ibot1;
  int ibot2;
  int itop;
  int itop1;
  int itop2;
  bool more;

  dtop = 0;
  dbot = 1;

  iarray = new int[n];
//
//  Compute the next permutation.
//
  more = false;

  for ( ; ; )
  {
    perm_next ( n, iarray, more, even );
//
//  The sign of this term depends on the sign of the permutation.
//
    if ( even )
    {
      itop = 1;
    }
    else
    {
      itop = -1;
    }
    ibot = 0;
//
//  Choose one item from each row, as specified by the permutation,
//  and multiply them.
//
    for ( i = 0; i < n; i++ )
    {
      itop1 = itop;
      ibot1 = ibot;
      itop2 = atop[i+(iarray[i]-1)*n];
      ibot2 = abot[i+(iarray[i]-1)*n];

      dec_mul ( itop1, ibot1, itop2, ibot2, dec_digit, itop, ibot );

    }
//
//  Add this term to the total.
//
    dec_add ( itop, ibot, dtop, dbot, dec_digit, dtop, dbot );

    if ( !more )
    {
      break;
    }
  }

  delete [] iarray;

  return;
}
//****************************************************************************80

void decmat_print ( int m, int n, int a[], int b[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    DECMAT_PRINT prints out decimal vectors and matrices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, int A[M*N], B[M*N], the decimal matrix.
//
//    Input, string TITLE, a title.
//
{
  int exponent;
  int i;
  int j;
  int jmax;
  int jmin;
  int kmax;
  int mantissa;
  int ncolum = 80;
  int npline;
  char *s;
//
//  Figure out how wide we must make each column.
//
  kmax = 0;

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      kmax = i4_max ( kmax, dec_width ( a[i+j*m], b[i+j*m] ) );
    }
  }

  s = new char[kmax];

  kmax = kmax + 2;
  npline = ncolum / kmax;
//
//  Now do the printing.
//
  for ( jmin = 0; jmin < n; jmin = jmin + npline )
  {
    jmax = i4_min ( jmin+npline-1, n-1 );

    cout << "\n";
    cout << title << "\n";
    cout << "\n";

    if ( 0 < jmin || jmax < n-1 )
    {
      cout << "Columns " << jmin << " to " << jmax << "\n";
      cout << "\n";
    }

    for ( i = 0; i < m; i++ )
    {
      for ( j = jmin; j <= jmax; j++ )
      {
        mantissa = a[i+j*m];
        exponent = b[i+j*m];
        s = dec_to_s ( mantissa, exponent );
        cout << setw(kmax) << s << "  ";
      }
      cout << "\n";
    }
  }

  delete [] s;

  return;
}
//****************************************************************************80

void derange_back_candidate ( int n, int a[], int k, int &nstack, int stack[], 
  int ncan[] )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGE_BACK_CANDIDATE finds possible values for the K-th entry of a derangement.
//
//  Discussion:
//
//    A derangement of N objects is a permutation which leaves no object
//    unchanged.
//
//    A derangement of N objects is a permutation with no fixed
//    points.  If we symbolize the permutation operation by "P",
//    then for a derangment, P(I) is never equal to I.
//
//    The number of derangements of N objects is sometimes called
//    the subfactorial function, or the derangement number D(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the derangement.
//
//    Input, int A[N].  The first K-1 entries of A
//    record the currently set values of the derangement.
//
//    Input, int K, the entry of the derangement for which candidates
//    are to be found.
//
//    Input/output, int &NSTACK, the length of the stack.
//
//    Input/output, int STACK[(N*(N+1))/2].  On output, we have added
//    the candidates for entry K to the end of the stack.
//
//    Input/output, int NCAN[N], the number of candidates for each level.
//
{
  int ican;
  int *ifree;
  int nfree;
//
//  Consider all the integers from 1 through N that have not been used yet.
//
  nfree = n - k + 1;
  ifree = new int[n];

  perm_free ( k-1, a, nfree, ifree );
//
//  Everything but K is a legitimate candidate for the K-th entry.
//
  ncan[k-1] = 0;

  for ( ican = 0; ican < nfree; ican++ )
  {
    if ( ifree[ican] != k )
    {
      ncan[k-1] = ncan[k-1] + 1;
      stack[nstack] = ifree[ican];
      nstack = nstack + 1;
    }

  }

  delete [] ifree;

  return;
}
//****************************************************************************80

void derange_back_next ( int n, int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGE_BACK_NEXT returns the next derangement of N items.
//
//  Discussion:
//
//    A derangement of N objects is a permutation which leaves no object
//    unchanged.
//
//    A derangement of N objects is a permutation with no fixed
//    points.  If we symbolize the permutation operation by "P",
//    then for a derangment, P(I) is never equal to I.
//
//    The number of derangements of N objects is sometimes called
//    the subfactorial function, or the derangement number D(N).
//
//    This routine uses backtracking.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items to be deranged.  N should be 2 or more.
//
//    Input/output, int A[N].
//    On the first call, the input value of A is not important.
//    On return with MORE = TRUE, A contains the next derangement.
//    On subsequent input, A should not be changed.
//
//    Input/output, bool &MORE.
//    On first call, set MORE to FALSE and do not alter it after.
//    On return, MORE is TRUE if another derangement is being returned in A,
//    and FALSE if no more derangements could be found.
//
{
  int i;
  static int indx = -1;
  static int k = -1;
  static int *ncan = NULL;
  static int *stack = NULL;
  static int stack_max = -1;
  static int stack_num = -1;

  if ( !( more ) )
  {
    if ( n < 2 )
    {
      more = false;
      return;
    }

    indx = 0;
    k = 0;
    stack_max = ( n * ( n + 1 ) ) / 2;
    stack_num = 0;

    if ( stack )
    {
      delete [] stack;
    }

    stack = new int[stack_max];
    for ( i = 0; i < stack_max; i++ )
    {
      stack[i] = 0;
    }

    if ( ncan )
    {
      delete [] ncan;
    }

    ncan = new int[n];
    for ( i = 0; i < n; i++ )
    {
      ncan[i] = 0;
    }
    more = true;
  }

  for ( ; ; )
  {
    i4vec_backtrack ( n, stack_max, stack, a, indx, k, stack_num, ncan );

    if ( indx == 1 )
    {
      break;
    }
    else if ( indx == 2 )
    {
      derange_back_candidate ( n, a, k, stack_num, stack, ncan );
    }
    else
    {
      more = false;
      delete [] ncan;
      ncan = NULL;
      delete [] stack;
      stack = NULL;
      break;
    }
  }

  return;
}
//****************************************************************************80

bool derange_check ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGE_CHECK is TRUE if a permutation is a derangement.
//
//  Discussion:
//
//    A derangement of N objects is a permutation which leaves no object
//    unchanged.
//
//    A derangement of N objects is a permutation with no fixed
//    points.  If we symbolize the permutation operation by "P",
//    then for a derangment, P(I) is never equal to I.
//
//    The number of derangements of N objects is sometimes called
//    the subfactorial function, or the derangement number D(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects permuted.
//
//    Input, int A[N], a permutation of the integers 1 through N.
//
//    Output, bool DERANGE_CHECK is TRUE if there was an error.
//
{
  int i;

  for ( i = 1; i <= n; i++ )
  {
    if ( a[i-1] == i )
    {
      return false;
    }
  }

  return true;
}
//****************************************************************************80

int derange_enum ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGE_ENUM returns the number of derangements of N objects.
//
//  Discussion:
//
//    A derangement of N objects is a permutation which leaves no object
//    unchanged.
//
//    A derangement of N objects is a permutation with no fixed
//    points.  If we symbolize the permutation operation by "P",
//    then for a derangment, P(I) is never equal to I.
//
//    The number of derangements of N objects is sometimes called
//    the subfactorial function, or the derangement number D(N).
//
//    D(N) is the number of ways of placing N non-attacking rooks on
//    an N by N chessboard with one diagonal deleted.
//
//    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
//
//    The number of permutations with exactly K items in the right
//    place is COMB(N,K) * D(N-K).
//
//    The formula:
//
//      D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
//
//    based on the inclusion/exclusion law.
//
//  Recursion:
//
//      D(0) = 1
//      D(1) = 0
//      D(2) = 1
//      D(N) = (N-1) * ( D(N-1) + D(N-2) )
//
//    or
//
//      D(0) = 1
//      D(1) = 0
//      D(N) = N * D(N-1) + (-1)**N
//
//  First values:
//
//     N         D(N)
//     0           1
//     1           0
//     2           1
//     3           2
//     4           9
//     5          44
//     6         265
//     7        1854
//     8       14833
//     9      133496
//    10     1334961
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Output, int DERANGE_ENUM, the number of derangements of N objects.
//
{
  int i;
  int value;
  int value1;
  int value2;

  if ( n < 0 )
  {
    value = 0;
  }
  else if ( n == 0 )
  {
    value = 1;
  }
  else if ( n == 1 )
  {
    value = 0;
  }
  else if ( n == 2 )
  {
    value = 1;
  }
  else
  {
    value1 = 0;
    value = 1;

    for ( i = 3; i <= n; i++ )
    {
      value2 = value1;
      value1 = value;
      value = ( i - 1 ) * ( value1 + value2 );
    }
  }

  return value;
}
//****************************************************************************80

void derange_enum2 ( int n, int d[] )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGE_ENUM2 returns the number of derangements of 0 through N objects.
//
//  Discussion:
//
//    A derangement of N objects is a permutation which leaves no object
//    unchanged.
//
//    A derangement of N objects is a permutation with no fixed
//    points.  If we symbolize the permutation operation by "P",
//    then for a derangment, P(I) is never equal to I.
//
//    The number of derangements of N objects is sometimes called
//    the subfactorial function, or the derangement number D(N).
//
//    D(N) is the number of ways of placing N non-attacking rooks on
//    an N by N chessboard with one diagonal deleted.
//
//    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
//
//    The number of permutations with exactly K items in the right
//    place is COMB(N,K) * D(N-K).
//
//    The formula is:
//
//      D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
//
//    based on the inclusion/exclusion law.
//
//  Recursion:
//
//      D(0) = 1
//      D(1) = 0
//      D(2) = 1
//      D(N) = (N-1) * ( D(N-1) + D(N-2) )
//
//    or
//
//      D(0) = 1
//      D(1) = 0
//      D(N) = N * D(N-1) + (-1)**N
//
//  Example:
//
//     N         D(N)
//     0           1
//     1           0
//     2           1
//     3           2
//     4           9
//     5          44
//     6         265
//     7        1854
//     8       14833
//     9      133496
//    10     1334961
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the maximum number of objects to be permuted.
//
//    Output, int D[N+1]; D(I) is the number of derangements of
//    I objects.
//
{
  int i;

  d[0] = 1;
  d[1] = 0;

  for ( i = 2; i <= n; i++ )
  {
    d[i] = ( i - 1 ) * ( d[i-1] + d[i-2] );
  }

  return;
}
//****************************************************************************80

int derange_enum3 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGE_ENUM3 returns the number of derangements of 0 through N objects.
//
//  Discussion:
//
//    A derangement of N objects is a permutation which leaves no object
//    unchanged.
//
//    A derangement of N objects is a permutation with no fixed
//    points.  If we symbolize the permutation operation by "P",
//    then for a derangment, P(I) is never equal to I.
//
//    The number of derangements of N objects is sometimes called
//    the subfactorial function, or the derangement number D(N).
//
//    D(N) is the number of ways of placing N non-attacking rooks on
//    an N by N chessboard with one diagonal deleted.
//
//    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
//
//    The number of permutations with exactly K items in the right
//    place is COMB(N,K) * D(N-K).
//
//    The formula is:
//
//      D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
//
//    based on the inclusion/exclusion law.
//
//    D(N) = nint ( N! / E )
//
//  Recursion:
//
//      D(0) = 1
//      D(1) = 0
//      D(2) = 1
//      D(N) = (N-1) * ( D(N-1) + D(N-2) )
//
//    or
//
//      D(0) = 1
//      D(1) = 0
//      D(N) = N * D(N-1) + (-1)**N
//
//  Example:
//
//     N         D(N)
//     0           1
//     1           0
//     2           1
//     3           2
//     4           9
//     5          44
//     6         265
//     7        1854
//     8       14833
//     9      133496
//    10     1334961
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the maximum number of objects to be permuted.
//
//    Output, int DERANGE_ENUM3, the number of derangements of N objects.
//
{
# define E 2.718281828459045

  int value;

  if ( n < 0 )
  {
    value = -1;
  }
  else if ( n == 0 )
  {
    value = 1;
  }
  else if ( n == 1 )
  {
    value = 0;
  }
  else
  {
    value = ( int ) ( 0.5 + ( r8_factorial ( n ) / E ) );
  }

  return value;
# undef E
}
//****************************************************************************80

void derange_weed_next ( int n, int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGE_WEED_NEXT computes all of the derangements of N objects, one at a time.
//
//  Discussion:
//
//    A derangement of N objects is a permutation which leaves no object
//    unchanged.
//
//    A derangement of N objects is a permutation with no fixed
//    points.  If we symbolize the permutation operation by "P",
//    then for a derangment, P(I) is never equal to I.
//
//    The number of derangements of N objects is sometimes called
//    the subfactorial function, or the derangement number D(N).
//
//    This routine simply generates all permutations, one at a time,
//    and weeds out those that are not derangements.
//
//  Example:
//
//    Here are the derangements when N = 4:
//
//    2143  3142  4123
//    2341  3412  4312
//    2413  3421  4321
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int A[N].
//    On first call, the input contents of A are unimportant.  But
//    on the second and later calls, the input value of A should be
//    the output value returned on the previous call.
//    On output, A contains the next derangement.
//
//    Input/output, bool &MORE.
//    Set MORE = FALSE before the first call.
//    MORE will be reset to TRUE and a derangement will be returned.
//    Each new call produces a new derangement until MORE is returned FALSE.
//
{
  bool deranged;
  static int maxder = 0;
  static int numder = 0;
//
//  Initialization on call with MORE = FALSE.
//
  if ( !( more ) )
  {
    maxder = derange_enum ( n );
    numder = 0;
  }
//
//  Watch out for cases where there are no derangements.
//
  if ( maxder == 0 )
  {
    more = false;
    return;
  }
//
//  Get the next permutation.
//
  for ( ; ; )
  {
    perm_lex_next ( n, a, more );
//
//  See if it is a derangment.
//
    deranged = derange_check ( n, a );

    if ( deranged )
    {
      break;
    }
  }

  numder = numder + 1;

  if ( maxder <= numder )
  {
    more = false;
  }

  return;
}
//****************************************************************************80

char digit_to_ch ( int digit )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_TO_CH returns the character representation of a decimal digit.
//
//  Example:
//
//    DIGIT   C
//    -----  ---
//      0    '0'
//      1    '1'
//    ...    ...
//      9    '9'
//     17    '*'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIGIT, the digit value between 0 and 9.
//
//    Output, char DIGIT_TO_CH, the corresponding character, or '*' if DIGIT
//    was illegal.
//
{
  if ( 0 <= digit && digit <= 9 )
  {
    return ( digit + 48 );
  }
  else
  {
    return '*';
  }
}
//****************************************************************************80

void digraph_arc_euler ( int nnode, int nedge, int inode[], int jnode[], 
  bool &success, int trail[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIGRAPH_ARC_EULER returns an Euler circuit in a digraph.
//
//  Discussion:
//
//    An Euler circuit of a digraph is a path which starts and ends at
//    the same node and uses each directed edge exactly once.  A digraph is
//    eulerian if it has an Euler circuit.  The problem is to decide whether
//    a given digraph is eulerian and to find an Euler circuit if the
//    answer is affirmative.
//
//  Method:
//
//    A digraph has an Euler circuit if and only if the number of incoming
//    edges is equal to the number of outgoing edges at each node.
//
//    This characterization gives a straightforward procedure to decide whether
//    a digraph is eulerian.  Furthermore, an Euler circuit in an eulerian
//    digraph G of NEDGE edges can be determined by the following method:
//
//      STEP 1: Choose any node U as the starting node, and traverse any edge
//      ( U, V ) incident to node U, and than traverse any unused edge incident
//      to node U.  Repeat this process of traversing unused edges until the
//      starting node U is reached.  Let P be the resulting walk consisting of
//      all used edges.  If all edges of G are in P, than stop.
//
//      STEP 2: Choose any unused edge ( X,  Y) in G such that X is
//      in P and Y is not in P.  Use node X as the starting node and
//      find another walk Q using all unused edges as in step 1.
//
//      STEP 3: Walk P and walk Q share a common node X, they can be merged
//      to form a walk R by starting at any node S of P and to traverse P
//      until node X is reached; than, detour and traverse all edges of Q
//      until node X is reached and continue to traverse the edges of P until
//      the starting node S is reached.  Set P = R.
//
//      STEP 4: Repeat steps 2 and 3 until all edges are used.
//
//    The running time of the algorithm is O ( NEDGE ).
//
//    The digraph is assumed to be connected.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Hang Tong Lau.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Hang Tong Lau,
//    Algorithms on Graphs,
//    Tab Books, 1989.
//
//  Parameters:
//
//    Input, int NNODE, the number of nodes.
//
//    Input, int NEDGE, the number of edges.
//
//    Input, int INODE[NEDGE], JNODE(NEDGE); the I-th edge starts at node
//    INODE(I) and ends at node JNODE(I).
//
//    Output, bool &SUCCESS, is TRUE if an Euler circuit was found.
//
//    Output, int TRAIL[NEDGE].  TRAIL[I] is the edge number of the I-th
//    edge in the Euler circuit.
//
{
  bool *candid;
  int *endnod;
  int i;
  int istak;
  int j;
  int k;
  int l;
  int len;
  int lensol;
  int lenstk;
  int *stack;
//
//  Check if the digraph is eulerian.
//
  for ( i = 0; i < nedge; i++ )
  {
    trail[i] = 0;
  }

  endnod = new int[nedge];

  for ( i = 0; i < nedge; i++ )
  {
    endnod[i] = 0;
  }

  for ( i = 1; i <= nedge; i++ )
  {
    j = inode[i-1];
    trail[j-1] = trail[j-1] + 1;
    j = jnode[i-1];
    endnod[j-1] = endnod[j-1] + 1;
  }

  for ( i = 1; i <= nnode; i++ )
  {
    if ( trail[i-1] != endnod[i-1] )
    {
      success = false;
      delete [] endnod;
      return;
    }
  }
//
//  The digraph is eulerian; find an Euler circuit.
//
  success = true;
  lensol = 1;
  lenstk = 0;

  candid = new bool[nedge];
  stack = new int[2*nedge];
//
//  Find the next edge.
//
  for ( ; ; )
  {
    if ( lensol == 1 )
    {
      endnod[0] = inode[0];
      stack[0] = 1;
      stack[1] = 1;
      lenstk = 2;
    }
    else
    {
      l = lensol - 1;

      if ( lensol != 2 )
      {
        endnod[l-1] = inode[trail[l-1]-1] + jnode[trail[l-1]-1] - endnod[l-2];
      }

      k = endnod[l-1];

      for ( i = 1; i <= nedge; i++ )
      {
        candid[i-1] = ( k == jnode[i-1] );
      }

      for ( i = 1; i <= l; i++ )
      {
        candid[trail[i-1]-1] = false;
      }

      len = lenstk;

      for ( i = 1; i <= nedge; i++ )
      {
        if ( candid[i-1] )
        {
          len = len + 1;
          stack[len-1] = i;
        }
      }
      stack[len] = len - lenstk;
      lenstk = len + 1;
    }

    for ( ; ; )
    {
      istak = stack[lenstk-1];
      lenstk = lenstk - 1;

      if ( istak != 0 )
      {
        break;
      }

      lensol = lensol - 1;

      if ( lensol == 0 )
      {
        i4vec_reverse ( nedge, trail );

        delete [] candid;
        delete [] endnod;
        delete [] stack;
        return;
      }
    }

    trail[lensol-1] = stack[lenstk-1];
    stack[lenstk-1] = istak - 1;

    if ( lensol == nedge )
    {
      break;
    }
    lensol = lensol + 1;
  }

  i4vec_reverse ( nedge, trail );

  delete [] candid;
  delete [] endnod;
  delete [] stack;

  return;
}
//****************************************************************************80

void digraph_arc_print ( int nedge, int inode[], int jnode[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    DIGRAPH_ARC_PRINT prints out a digraph from an edge list.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NEDGE, the number of edges.
//
//    Input, int INODE[NEDGE], JNODE[NEDGE], the beginning and end
//    nodes of the edges.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  for ( i = 0; i < nedge; i++ )
  {
    cout << setw(6) << i+1 << "  "
         << setw(6) << inode[i] << "  "
         << setw(6) << jnode[i] << "\n";
  }

  return;
}
//****************************************************************************80

void diophantine ( int a, int b, int c, bool &error, int &x, int &y )

//****************************************************************************80
//
//  Purpose:
//
//    DIOPHANTINE solves a Diophantine equation A * X + B * Y = C.
//
//  Discussion:
//
//    Given integers A, B and C, produce X and Y so that
//
//      A * X + B * Y = C.
//
//    In general, the equation is solvable if and only if the
//    greatest common divisor of A and B also divides C.
//
//    A solution (X,Y) of the Diophantine equation also gives the solution
//    X to the congruence equation:
//
//      A * X = C mod ( B ).
//
//    Generally, if there is one nontrivial solution, there are an infinite
//    number of solutions to a Diophantine problem.
//    If (X0,Y0) is a solution, then so is ( X0+T*B/D, Y0-T*A/D ) where
//    T is any integer, and D is the greatest common divisor of A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric Weisstein, editor,
//    CRC Concise Encylopedia of Mathematics,
//    CRC Press, 1998, page 446.
//
//  Parameters:
//
//    Input, int A, B, C, the coefficients of the Diophantine equation.
//
//    Output, bool &ERROR, is TRUE if an error occurred.
//
//    Output, int &X, &Y, the solution of the Diophantine equation.
//    Note that the algorithm will attempt to return a solution with
//    smallest Euclidean norm.
//
{
# define N_MAX 100

  int a_copy;
  int a_mag;
  int a_sign;
  int b_copy;
  int b_mag;
  int b_sign;
  int c_copy;
  int g;
  int k;
  int n;
  double norm_new;
  double norm_old;
  int q[N_MAX];
  bool swap;
  int temp;
  int xnew;
  int ynew;
//
//  Defaults for output parameters.
//
  error = false;
  x = 0;
  y = 0;
//
//  Special cases.
//
  if ( a == 0 && b == 0 && c == 0 )
  {
    x = 0;
    y = 0;
    return;
  }
  else if ( a == 0 && b == 0 && c != 0 )
  {
    error = true;
    x = 0;
    y = 0;
    return;
  }
  else if ( a == 0 && b != 0 && c == 0 )
  {
    x = 0;
    y = 0;
    return;
  }
  else if ( a == 0 && b != 0 && c != 0 )
  {
    x = 0;
    y = c / b;
    if ( ( c % b ) != 0 )
    {
      error = true;
    }
    return;
  }
  else if ( a != 0 && b == 0 && c == 0 )
  {
    x = 0;
    y = 0;
    return;
  }
  else if ( a != 0 && b == 0 && c != 0 )
  {
    x = c / a;
    y = 0;
    if ( ( c % a ) != 0 )
    {
      error = true;
    }
    return;
  }
  else if ( a != 0 && b != 0 && c == 0 )
  {
    g = i4_gcd ( a, b );
    x = b / g;
    y = - a / g;
    return;
  }
//
//  Now handle the "general" case: A, B and C are nonzero.
//
//  Step 1: Compute the GCD of A and B, which must also divide C.
//
  g = i4_gcd ( a, b );

  if ( ( c % g ) != 0 )
  {
    error = true;
    return;
  }

  a_copy = a / g;
  b_copy = b / g;
  c_copy = c / g;
//
//  Step 2: Split A and B into sign and magnitude.
//
  a_mag = abs ( a_copy );
  a_sign = i4_sign ( a_copy );
  b_mag = abs ( b_copy );
  b_sign = i4_sign ( b_copy );
//
//  Another special case, A_MAG = 1 or B_MAG = 1.
//
  if ( a_mag == 1 )
  {
    x = a_sign * c_copy;
    y = 0;
    return;
  }
  else if ( b_mag == 1 )
  {
    x = 0;
    y = b_sign * c_copy;
    return;
  }
//
//  Step 3: Produce the Euclidean remainder sequence.
//
  if ( b_mag <= a_mag )
  {
    swap = false;
    q[0] = a_mag;
    q[1] = b_mag;
  }
  else
  {
    swap = true;
    q[0] = b_mag;
    q[1] = a_mag;
  }

  n = 3;

  for ( ; ; )
  {
    q[n-1] = q[n-3] % q[n-2];

    if ( q[n-1] == 1 )
    {
      break;
    }

    n = n + 1;

    if ( N_MAX < n )
    {
      error = true;
      cerr << "\n";
      cerr << "DIOPHANTINE - Fatal error!\n";
      cerr << "  Exceeded number of iterations.\n";
      exit ( 1 );
    }
  }
//
//  Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
//
  y = 0;
  for ( k = n; 2 <= k; k-- )
  {
    x = y;
    y = ( 1 - x * q[k-2] ) / q[k-1];
  }
//
//  Step 5: Undo the swapping.
//
  if ( swap )
  {
    i4_swap ( x, y );
  }
//
//  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
//
  x = x * a_sign;
  y = y * b_sign;
//
//  Step 7: Multiply by C, so that X * A + Y * B = C.
//
  x = x * c_copy;
  y = y * c_copy;
//
//  Step 8: Given a solution (X,Y), try to find the solution of
//  minimal magnitude.
//
  diophantine_solution_minimize ( a_copy, b_copy, x, y );

  return;
# undef N_MAX
}
//****************************************************************************80

void diophantine_solution_minimize ( int a, int b, int &x, int &y )

//****************************************************************************80
//
//  Purpose:
//
//    DIOPHANTINE_SOLUTION_MINIMIZE seeks a minimal solution of a Diophantine equation.
//
//  Discussion:
//
//    Given a solution (X,Y) of a Diophantine equation:
//
//      A * X + B * Y = C.
//
//    then there are an infinite family of solutions of the form
//
//      ( X(i), Y(i) ) = ( X + i * B, Y - i * A )
//
//    An integral solution of minimal Euclidean norm can be found by
//    tentatively moving along the vectors (B,-A) and (-B,A) one step
//    at a time.
//
//    When large integer values are input, the real arithmetic used
//    is essential.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric Weisstein, editor,
//    CRC Concise Encylopedia of Mathematics,
//    CRC Press, 1998, page 446.
//
//  Parameters:
//
//    Input, int A, B, the coefficients of the Diophantine equation.
//    A and B are assumed to be relatively prime.
//
//    Input/output, int &X, &Y, on input, a solution of the Diophantine
//    equation.  On output, a solution of minimal Euclidean norm.
//
{
  double fa;
  double fb;
  double fx;
  double fy;
  double norm;
  double norm_new;
  double t;
  int xnew;
  int ynew;
//
//  Compute the minimum for T real, and then look nearby.
//
  fa = ( double ) a;
  fb = ( double ) b;
  fx = ( double ) x;
  fy = ( double ) y;

  t = ( - fb * fx + fa * fy ) / ( fa * fa + fb * fb );

  x = x + r8_nint ( t ) * b;
  y = y - r8_nint ( t ) * a;
//
//  Now look nearby.
//
  norm = ( fx * fx + fy * fy );

  for ( ; ; )
  {
    xnew = x + b;
    ynew = y - a;

    fx = ( double ) xnew;
    fy = ( double ) ynew;

    norm_new = ( fx * fx + fy * fy );

    if ( norm <= norm_new )
    {
      break;
    }

    x = xnew;
    y = ynew;
    norm = norm_new;
  }

  for ( ; ; )
  {
    xnew = x - b;
    ynew = y + a;

    fx = ( double ) xnew;
    fy = ( double ) ynew;

    norm_new = ( fx * fx + fy * fy );

    if ( norm <= norm_new )
    {
      break;
    }

    x = xnew;
    y = ynew;
    norm = norm_new;
  }

  return;
}
//****************************************************************************80

void dvec_add ( int n, int dvec1[], int dvec2[], int dvec3[] )

//****************************************************************************80
//
//  Purpose:
//
//    DVEC_ADD adds two (signed) decimal vectors.
//
//  Discussion:
//
//    A DVEC is an integer vector of decimal digits, intended to
//    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
//    is the coefficient of 10**(N-2), and DVEC(N) contains sign
//    information.  It is 0 if the number is positive, and 9 if
//    the number is negative.
//
//  Example:
//
//    N = 4
//
//      DVEC1     +   DVEC2     =   DVEC3
//
//    ( 0 0 1 7 ) + ( 0 1 0 4 ) = ( 0 0 1 2 1 )
//
//          17    +       104   =         121
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int DVEC1[N], DVEC2[N], the vectors to be added.
//
//    Output, int DVEC3[N], the sum of the two input vectors.
//
{
  int base = 10;
  int i;
  bool overflow;

  overflow = false;

  for ( i = 0; i < n; i++ )
  {
    dvec3[i] = dvec1[i] + dvec2[i];
  }

  for ( i = 0; i < n; i++ )
  {
    while ( base <= dvec3[i] )
    {
      dvec3[i] = dvec3[i] - base;
      if ( i < n-1 )
      {
        dvec3[i+1] = dvec3[i+1] + 1;
      }
      else
      {
        overflow = true;
        return;
      }
    }
  }

  return;
}
//****************************************************************************80

void dvec_complementx ( int n, int dvec1[], int dvec2[] )

//****************************************************************************80
//
//  Purpose:
//
//    DVEC_COMPLEMENTX computes the ten's complement of a decimal vector.
//
//  Discussion:
//
//    A DVEC is an integer vector of decimal digits, intended to
//    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
//    is the coefficient of 10**(N-2), and DVEC(N) contains sign
//    information.  It is 0 if the number is positive, and 9 if
//    the number is negative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int DVEC1[N], the vector to be complemented.
//
//    Output, int DVEC2[N], the complemented vector.
//
{
  int base = 10;
  int *dvec3;
  int *dvec4;
  int i;

  dvec3 = new int[n];
  dvec4 = new int[n];
  
  for ( i = 0; i < n; i++ )
  {
    dvec3[i] = ( base - 1 ) - dvec1[i];
  }

  dvec4[0] = 1;
  for ( i = 1; i < n; i++ )
  {
    dvec4[i] = 0;
  }

  dvec_add ( n, dvec3, dvec4, dvec2 );

  delete [] dvec3;
  delete [] dvec4;

  return;
}
//****************************************************************************80

void dvec_mul ( int n, int dvec1[], int dvec2[], int dvec3[] )

//****************************************************************************80
//
//  Purpose:
//
//    DVEC_MUL computes the product of two decimal vectors.
//
//  Discussion:
//
//    A DVEC is an integer vector of decimal digits, intended to
//    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
//    is the coefficient of 10**(N-2), and DVEC(N) contains sign
//    information.  It is 0 if the number is positive, and 9 if
//    the number is negative.
//
//    Since the user may want to make calls like
//
//      dvec_mul ( n, dvec1, dvec1, dvec3 )
//    or even
//      dvec_mul ( n, dvec1, dvec1, dvec1 )
//
//    we need to copy the arguments, work on them, and then copy out the result.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int DVEC1[N], DVEC2[N], the vectors to be multiplied.
//
//    Output, int DVEC3[N], the product of the two input vectors.
//
{
  int base = 10;
  int carry;
  int *dveca;
  int *dvecb;
  int *dvecc;
  int i;
  int j;
  int product_sign;

  dveca = new int[n];
  dvecb = new int[n];
  dvecc = new int[n];
//
//  Copy the input.
//
  for ( i = 0; i < n; i++ )
  {
    dveca[i] = dvec1[i];
  }

  for ( i = 0; i < n; i++ )
  {
    dvecb[i] = dvec2[i];
  }
//
//  Record the sign of the product.
//  Make the factors positive.
//
  product_sign = 1;

  if ( dveca[n-1] != 0 )
  {
    product_sign = - product_sign;
    dvec_complementx ( n, dveca, dveca );
  }

  if ( dvecb[n-1] != 0 )
  {
    product_sign = - product_sign;
    dvec_complementx ( n, dvecb, dvecb );
  }

  for ( i = 0; i < n; i++ )
  {
    dvecc[i] = 0;
  }
//
//  Multiply.
//
  for ( i = 1; i <= n-1; i++ )
  {
    for ( j = i; j <= n-1; j++ )
    {
      dvecc[j-1] = dvecc[j-1] + dveca[i-1] * dvecb[j-i];
    }
  }
//
//  Take care of carries.
//  Unlike the DVEC_ADD routine, we do NOT allow carries into the
//  N-th position.
//
  for ( i = 0; i < n-1; i++ )
  {
    carry = dvecc[i] / base;
    dvecc[i] = dvecc[i] - carry * base;

    if ( i < n - 2 )
    {
      dvecc[i+1] = dvecc[i+1] + carry;
    }
  }
//
//  Take care of the sign of the product.
//
  if ( product_sign < 0 )
  {
    dvec_complementx ( n, dvecc, dvecc );
  }
//
//  Copy the output.
//
  for ( i = 0; i < n; i++ )
  {
    dvec3[i] = dvecc[i];
  }

  delete [] dveca;
  delete [] dvecb;
  delete [] dvecc;

  return;
}
//****************************************************************************80

void dvec_print ( int n, int dvec[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    DVEC_PRINT prints a decimal integer vector, with an optional title.
//
//  Discussion:
//
//    A DVEC is an integer vector of decimal digits, intended to
//    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
//    is the coefficient of 10**(N-2), and DVEC(N) contains sign
//    information.  It is 0 if the number is positive, and 9 if
//    the number is negative.
//
//    The vector is printed "backwards", that is, the first entry
//    printed is DVEC(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int DVEC[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int ihi;
  int ilo;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  for ( ihi = n; 1 <= ihi; ihi = ihi - 80 )
  {
    cout << "  ";
    ilo = i4_max ( ihi - 80 + 1, 1 );
    for ( i = ihi; ilo <= i; i-- )
    {
      cout << dvec[i-1];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void dvec_sub ( int n, int dvec1[], int dvec2[], int dvec3[] )

//****************************************************************************80
//
//  Purpose:
//
//    DVEC_SUB subtracts two decimal vectors.
//
//  Discussion:
//
//    A DVEC is an integer vector of decimal digits, intended to
//    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
//    is the coefficient of 10**(N-2), and DVEC(N) contains sign
//    information.  It is 0 if the number is positive, and 9 if
//    the number is negative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the length of the vectors.
//
//    Input, int DVEC1[N], DVEC2[N]), the vectors to be subtracted.
//
//    Output, int DVEC3[N], the value of DVEC1 - DVEC2.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    dvec3[i] = dvec2[i];
  }

  dvec_complementx ( n, dvec3, dvec3 );

  dvec_add ( n, dvec1, dvec3, dvec3 );

  return;
}
//****************************************************************************80

int dvec_to_i4 ( int n, int dvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    DVEC_TO_I4 makes an integer from a (signed) decimal vector.
//
//  Discussion:
//
//    A DVEC is an integer vector of decimal digits, intended to
//    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
//    is the coefficient of 10**(N-2), and DVEC(N) contains sign
//    information.  It is 0 if the number is positive, and 9 if
//    the number is negative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vector.
//
//    Input, int DVEC[N], the decimal vector.
//
//    Output, int DVEC_TO_I4, the integer.
//
{
  int base = 10;
  int *dvec2;
  int i;
  int i_sign;
  int i4;

  dvec2 = new int[n];

  for ( i = 0; i < n; i++ )
  {
    dvec2[i] = dvec[i];
  }

  i_sign = 1;

  if ( dvec2[n-1] == base - 1 )
  {
    i_sign = -1;
    dvec_complementx ( n-1, dvec2, dvec2 );
  }

  i4 = 0;
  for ( i = n-2; 0 <= i4; i4-- )
  {
    i4 = base * i4 + dvec2[i];
  }

  i4 = i_sign * i4;

  delete [] dvec2;

  return i4;
}
//****************************************************************************80

void equiv_next ( int n, int &npart, int jarray[], int iarray[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    EQUIV_NEXT computes the partitions of a set one at a time.
//
//  Discussion:
//
//    A partition of a set assigns each element to exactly one subset.
//
//    The number of partitions of a set of size N is the Bell number B(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, number of elements in the set to be partitioned.
//
//    Output, int &NPART, number of subsets in the partition.
//
//    Output, int JARRAY[N].  JARRAY[I] is the number of elements
//    in the I-th subset of the partition.
//
//    Output, int IARRAY[N].  IARRAY(I) is the class to which
//    element I belongs.
//
//    Input/output, bool &MORE.  Set MORE = FALSE before first call.
//    It is reset and held at TRUE as long as
//    the partition returned is not the last one.
//    When MORE is returned FALSE, all the partitions
//    have been computed and returned.
//
{
  int i;
  int l;
  int m;

  if ( !more )
  {
    npart = 1;
    for ( i = 0; i < n; i++ )
    {
      iarray[i] = 1;
    }
    jarray[0] = n;
  }
  else
  {
    m = n;

    while ( jarray[iarray[m-1]-1] == 1 )
    {
      iarray[m-1] = 1;
      m = m - 1;
    }

    l = iarray[m-1];
    npart = npart + m - n;
    jarray[0] = jarray[0] + n - m;

    if ( l == npart )
    {
      npart = npart + 1;
      jarray[npart-1] = 0;
    }
    iarray[m-1] = l + 1;
    jarray[l-1] = jarray[l-1] - 1;
    jarray[l] = jarray[l] + 1;
  }

  more = ( npart != n );

  return;
}
//****************************************************************************80

void equiv_next2 ( bool &done, int iarray[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    EQUIV_NEXT2 computes, one at a time, the partitions of a set.
//
//  Discussion:
//
//    A partition of a set assigns each element to exactly one subset.
//
//    The number of partitions of a set of size N is the Bell number B(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2003
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input/output, bool &DONE.  Before the very first call, the
//    user should set DONE to TRUE, which prompts the program
//    to initialize its data, and return the first partition.
//    Thereafter, the user should call again, for the next
//    partition, and so on, until the routine returns with DONE
//    equal to TRUE, at which point there are no more partitions
//    to compute.
//
//    Input/output, int IARRAY[N], contains the information
//    defining the current partition.  The user should not alter
//    IARRAY between calls.  Except for the very first
//    call, the routine uses the previous output value of IARRAY to compute
//    the next value.
//    The entries of IARRAY are the partition subset to which each
//    element of the original set belongs.  If there are NPART distinct
//    parts of the partition, then each entry of IARRAY will be a
//    number between 1 and NPART.  Every number from 1 to NPART will
//    occur somewhere in the list.  If the entries of IARRAY are
//    examined in order, then each time a new partition subset occurs,
//    it will be the next unused integer.
//    For instance, for N = 4, the program will describe the set
//    where each element is in a separate subset as 1, 2, 3, 4,
//    even though such a partition might also be described as
//    4, 3, 2, 1 or even 1, 5, 8, 19.
//
//    Input, int N, the number of elements in the set.
//
{
  int i;
  int imax;
  int j;
  int jmax;

  if ( done )
  {
    done = false;
    for ( i = 0; i < n; i++ )
    {
      iarray[i] = 1;
    }
  }
  else
  {
//
//  Find the last element J that can be increased by 1.
//  This is the element that is not equal to its maximum possible value,
//  which is the maximum value of all preceding elements +1.
//
    jmax = iarray[0];
    imax = 1;

    for ( j = 2; j <= n; j++ )
    {
      if ( jmax < iarray[j-1] )
      {
        jmax = iarray[j-1];
      }
      else
      {
        imax = j;
      }
    }
//
//  If no element can be increased by 1, we are done.
//
    if ( imax == 1 )
    {
      done = true;
      return;
    }
//
//  Increase the value of the IMAX-th element by 1, set its successors to 1.
//
    done = false;
    iarray[imax-1] = iarray[imax-1] + 1;
    for ( i = imax; i < n; i++ )
    {
      iarray[i] = 1;
    }
  }
  return;
}
//****************************************************************************80

void equiv_print ( int n, int iarray[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    EQUIV_PRINT prints a partition of a set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, number of elements in set to be partitioned.
//
//    Input, int IARRAY[N], defines the partition or set of equivalence
//    classes.  Element I belongs to subset IARRAY[I].
//
//    Input, string TITLE, a title.
//
{
  int *karray;
  int j;
  int k;
  int kk;
  int s;
  int s_max;
  int s_min;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  cout << "  Set    Size  Elements\n";

  karray = new int[n];

  s_min = i4vec_min ( n, iarray );
  s_max = i4vec_max ( n, iarray );

  for ( s = s_min; s <= s_max; s++ )
  {
    k = 0;

    for ( j = 0; j < n; j++ )
    {
      if ( iarray[j] == s )
      {
        karray[k] = j+1;
        k = k + 1;
      }
    }

    if ( 0 < k )
    {
      cout                 << "  "
           << setw(4) << s << "  "
           << setw(4) << k << " :: ";
      for ( kk = 0; kk < k; kk++ )
      {
        cout << setw(4) << karray[kk] << "  ";
      }
      cout << "\n";
    }
  }

  delete [] karray;

  return;
}
//****************************************************************************80

void equiv_print2 ( int n, int s[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    EQUIV_PRINT2 prints a partition of a set.
//
//  Discussion:
//
//    The partition is printed using the parenthesis format.
//
//    For example, here are the partitions of a set of 4 elements:
//
// ÊÊÊÊÊ(1,2,3,4)
// ÊÊÊÊÊ(1,2,3)(4)
// ÊÊÊÊÊ(1,2,4)(3)
// ÊÊÊÊÊ(1,2)(3,4)
// ÊÊÊÊÊ(1,2)(3)(4)
// ÊÊÊÊÊ(1,3,4)(2)
// ÊÊÊÊÊ(1,3)(2,4)
// ÊÊÊÊÊ(1,3)(2)(4)
// ÊÊÊÊÊ(1,4)(2,3)
// ÊÊÊÊÊ(1)(2,3,4)
// ÊÊÊÊÊ(1)(2,3)(4)
// ÊÊÊÊÊ(1,4)(2)(3)
// ÊÊÊÊÊ(1)(2,4)(3)
// ÊÊÊÊÊ(1)(2)(3,4)
// ÊÊÊÊÊ(1)(2)(3)(4)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, number of elements in the set.
//
//    Input, int S[N], defines the partition.  
//    Element I belongs to subset S[I].
//
//    Input, string TITLE, a title.
//
{
  int i;
  int j;
  int s_max;
  int s_min;
  int size;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  s_min = i4vec_min ( n, s );
  s_max = i4vec_max ( n, s );

  for ( j = s_min; j <= s_max; j++ )
  {
    cout << "(";
    size = 0;
    for ( i = 0; i < n; i++ )
    {
      if ( s[i] == j )
      {
        if ( 0 < size )
        {
          cout << ",";
        }
        cout << i;
        size = size + 1;
      }
    }
    cout << ")";
  }
  cout << "\n";

  return;
}
//****************************************************************************80

void equiv_random ( int n, int &seed, int &npart, int a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    EQUIV_RANDOM selects a random partition of a set.
//
//  Discussion:
//
//    The user does not control the number of parts in the partition.
//
//    The equivalence classes are numbered in no particular order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of elements in the set to be partitioned.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int &NPART, the number of classes or parts in the 
//    partition.  NPART will be between 1 and N.
//
//    Output, int A[N], indicates the class to which each element
//    is assigned.
//
//    Output, double B[N].  B(K) = C(K)/(K!), where
//    C(K) = number of partitions of a set of K objects.
//
{
  int k;
  int l;
  int m;
  double sum1;
  double z;
  double zhi = 1.0;
  double zlo = 0.0;

  b[0] = 1.0;

  for ( l = 1; l <= n-1; l++ )
  {
    sum1 = 1.0 / ( double ) l;
    for ( k = 1; k <= l-1; k++ )
    {
      sum1 = ( sum1 + b[k-1] ) / ( double ) ( l - k );
    }

    b[l] = ( sum1 + b[l-1] ) / ( double ) ( l + 1 );
  }

  m = n;
  npart = 0;

  for ( ; ; )
  {
    z = r8_uniform ( zlo, zhi, seed );
    z = ( double ) ( m ) * b[m-1] * z;
    k = 0;
    npart = npart + 1;

    while ( 0.0 <= z )
    {
      a[m-1] = npart;
      m = m - 1;

      if ( m == 0 )
      {
        break;
      }

      z = z - b[m-1];
      k = k + 1;
      z = z * k;
    }

    if ( m == 0 )
    {
      break;
    }
  }
//
//  Randomly permute the assignments.
//
  perm_random2 ( n, seed, a );

  return;
}
//****************************************************************************80

void euler ( int n, int ieuler[] )

//****************************************************************************80
//
//  Purpose:
//
//    EULER returns the N-th row of Euler's triangle.
//
//  Discussion:
//
//    E(N,K) counts the number of permutations of the N digits that have
//    exactly K "ascents", that is, K places where the Ith digit is
//    less than the (I+1)th digit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the row of Euler's triangle desired.
//
//    Output, int IEULER[N+1], the N-th row of Euler's
//    triangle, IEULER[K] contains the value of E(N,K).  Note
//    that IEULER[0] should be 1 and IEULER[N] should be 0.
//
{
  int irow;
  int k;

  ieuler[0] = 1;

  if ( 0 < n )
  {
    ieuler[1] = 0;

    for ( irow = 2; irow <= n; irow++ )
    {
      ieuler[irow] = 0;

      for ( k = irow-1; 1 <= k; k-- )
      {
        ieuler[k] = ( k + 1 ) * ieuler[k] + ( irow - k ) * ieuler[k-1];
      }
      ieuler[0] = 1;
    }
  }
  return;
}
//****************************************************************************80

int frobenius_number_order2 ( int c1, int c2 )

//****************************************************************************80
//
//  Purpose:
//
//    FROBENIUS_NUMBER_ORDER2 returns the Frobenius number for order 2.
//
//  Discussion:
//
//    The Frobenius number of order N is the solution of the Frobenius
//    coin sum problem for N coin denominations.
//
//    The Frobenius coin sum problem assumes the existence of
//    N coin denominations, and asks for the largest value that cannot
//    be formed by any combination of coins of these denominations.
//
//    The coin denominations are assumed to be distinct positive integers.
//
//    For general N, this problem is fairly difficult to handle.
//
//    For N = 2, it is known that:
//
//    * if C1 and C2 are not relatively prime, then
//      there are infinitely large values that cannot be formed.
//
//    * otherwise, the largest value that cannot be formed is
//      C1 * C2 - C1 - C2, and that exactly half the values between
//      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
//
//    As a simple example, if C1 = 2 and C2 = 7, then the largest
//    unrepresentable value is 5, and there are (5+1)/2 = 3
//    unrepresentable values, namely 1, 3, and 5.
//
//    For a general N, and a set of coin denominations C1, C2, ..., CN,
//    the Frobenius number F(N, C(1:N) ) is defined as the largest value
//    B for which the equation
//
//      C1*X1 + C2*X2 + ... + CN*XN = B
//
//    has no nonnegative integer solution X(1:N).
//
//    In the Mathematica Package "NumberTheory", the Frobenius number
//    can be determined by
//
//    <<NumberTheory`Frobenius`
//    FrobeniusF[ {C1,...,CN} ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Sylvester,
//    Question 7382,
//    Mathematical Questions with their Solutions,
//    Educational Times,
//    Volume 41, page 21, 1884.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input, int C1, C2, the coin denominations. C1 and C2
//    should be positive and relatively prime.
//
//    Output, int FROBENIUS_NUMBER_ORDER2, the Frobenius number of (C1,C2).
//
{
  int value;

  if ( c1 <= 0 )
  {
    value = i4_huge ( );
  }
  else if ( c2 <= 0 )
  {
    value = i4_huge ( );
  }
  else if ( i4_gcd ( c1, c2 ) != 1 )
  {
    value = i4_huge ( );
  }
  else
  {
    value = c1 * c2 - c1 - c2;
  }

  return value;
}
//****************************************************************************80

void frobenius_number_order2_values ( int &n_data, int &c1, int &c2, int &f )

//****************************************************************************80
//
//  Purpose:
//
//    FROBENIUS_NUMBER_ORDER2_VALUES returns values of the order 2 Frobenius number.
//
//  Discussion:
//
//    The Frobenius number of order N is the solution of the Frobenius
//    coin sum problem for N coin denominations.
//
//    The Frobenius coin sum problem assumes the existence of
//    N coin denominations, and asks for the largest value that cannot
//    be formed by any combination of coins of these denominations.
//
//    The coin denominations are assumed to be distinct positive integers.
//
//    For general N, this problem is fairly difficult to handle.
//
//    For N = 2, it is known that:
//
//    * if C1 and C2 are not relatively prime, then
//      there are infinitely large values that cannot be formed.
//
//    * otherwise, the largest value that cannot be formed is
//      C1 * C2 - C1 - C2, and that exactly half the values between
//      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
//
//    As a simple example, if C1 = 2 and C2 = 7, then the largest
//    unrepresentable value is 5, and there are (5+1)/2 = 3
//    unrepresentable values, namely 1, 3, and 5.
//
//    For a general N, and a set of coin denominations C1, C2, ..., CN,
//    the Frobenius number F(N, C(1:N) ) is defined as the largest value
//    B for which the equation
//
//      C1*X1 + C2*X2 + ... + CN*XN = B
//
//    has no nonnegative integer solution X(1:N).
//
//    In the Mathematica Package "NumberTheory", the Frobenius number
//    can be determined by
//
//    <<NumberTheory`Frobenius`
//    FrobeniusF[ {C1,...,CN} ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Sylvester,
//    Question 7382,
//    Mathematical Questions with their Solutions,
//    Educational Times,
//    Volume 41, page 21, 1884.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &C1, &C2, the parameters of the function.
//
//    Output, int &F, the value of the function.
//
{
# define N_MAX 6

  int c1_vec[N_MAX] = {
     2, 
     3, 
     4, 
     5, 
    12, 
    99 };
  int c2_vec[N_MAX] = {
      5, 
     17, 
     19, 
     13, 
     11, 
    100 };
  int f_vec[N_MAX] = {
    3, 
    31, 
    23, 
    47, 
    109, 
    9701 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    c1 = 0;
    c2 = 0;
    f = 0;
  }
  else
  {
    c1 = c1_vec[n_data-1];
    c2 = c2_vec[n_data-1];
    f = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gamma_log_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_LOG_VALUES returns some values of the Log Gamma function for testing.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int &N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA isincremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 18

  double bvec[N_MAX] = {
     1.524064183E+00,    0.7966780066E+00,   0.3982337117E+00,  
     0.1520599127E+00,   0.000000000E+00,   -0.04987246543E+00, 
    -0.08537410945E+00, -0.1081747934E+00,  -0.1196128950E+00,  
    -0.1207822040E+00,  -0.1125917658E+00,  -0.09580771625E+00, 
    -0.07108385116E+00, -0.03898428380E+00,  0.000000000E+00,   
    12.80182743E+00,    39.33988571E+00,    71.25704193E+00 };
  double xvec[N_MAX] = {
    0.2E+00,  0.4E+00,  0.6E+00,  0.8E+00, 
    1.0E+00,  1.1E+00,  1.2E+00,  1.3E+00, 
    1.4E+00,  1.5E+00,  1.6E+00,  1.7E+00, 
    1.8E+00,  1.9E+00,  2.0E+00, 10.0E+00, 
   20.0E+00, 30.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  if ( N_MAX <= n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = xvec[n_data];
    fx = bvec[n_data];
    n_data = n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

unsigned long get_seed ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, unsigned long GET_SEED, a random seed value.
//
{
# define ULONG_MAX2 4294967295UL

  time_t clock;
  int i;
  int hours;
  int minutes;
  int seconds;
  struct tm *lt;
  static unsigned long seed = 0;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  if ( seed == 0 )
  {
    clock = time ( &tloc );
    lt = localtime ( &clock );
//
//  Extract HOURS.
//
    hours = lt->tm_hour;
//
//  In case of 24 hour clocks, shift so that HOURS is between 1 and 12.
//
    if ( 12 < hours )
    {
      hours = hours - 12;
    }
//
//  Move HOURS to 0, 1, ..., 11
//
    hours = hours - 1;

    minutes = lt->tm_min;

    seconds = lt->tm_sec;

    seed = seconds + 60 * ( minutes + 60 * hours );
//
//  We want values in [1,43200], not [0,43199].
//
    seed = seed + 1;
//
//  Remap SEED from [1,43200] to [1,ULONG_MAX2].
//
    seed = ( unsigned long ) 
      ( ( ( double ) seed )
      * ( ( double ) ULONG_MAX2 ) / ( 60.0 * 60.0 * 12.0 ) );
  }
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;

# undef ULONG_MAX2
}
//****************************************************************************80

void gray_next ( int n, int &change )

//****************************************************************************80
//
//  Purpose:
//
//    GRAY_NEXT generates the next Gray code by switching one item at a time.
//
//  Discussion:
//
//    On the first call only, the user must set CHANGE = -N.
//    This initializes the routine to the Gray code for N zeroes.
//
//    Each time it is called thereafter, it returns in CHANGE the index
//    of the item to be switched in the Gray code.  The sign of CHANGE
//    indicates whether the item is to be added or subtracted (or
//    whether the corresponding bit should become 1 or 0).  When
//    CHANGE is equal to N+1 on output, all the Gray codes have been
//    generated.
//
//    The routine has internal memory that is set up on call with
//    CHANGE = -N, and released on final return.
//
//  Example:
//
//    N  CHANGE         Subset in/out   Binary Number
//                      Interpretation  Interpretation
//                       1 2 4 8
//   --  ---------      --------------  --------------
//
//    4   -4 / 0         0 0 0 0         0
//
//        +1             1 0 0 0         1
//          +2           1 1 0 0         3
//        -1             0 1 0 0         2
//            +3         0 1 1 0         6
//        +1             1 1 1 0         7
//          -2           1 0 1 0         5
//        -1             0 0 1 0         4
//              +4       0 0 1 1        12
//        +1             1 0 1 1        13
//          +2           1 1 1 1        15
//        -1             0 1 1 1        14
//            -3         0 1 0 1        10
//        +1             1 1 0 1        11
//          -2           1 0 0 1         9
//        -1             0 0 0 1         8
//              -4       0 0 0 0         0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the order of the total set from which
//    subsets will be drawn.
//
//    Input/output, int &CHANGE.  This item is used for input only
//    on the first call for a particular sequence of Gray codes,
//    at which time it must be set to -N.  This corresponds to
//    all items being excluded, or all bits being 0, in the Gray code.
//    On output, CHANGE indicates which of the N items must be "changed", 
//    and the sign indicates whether the item is to be added or removed
//    (or the bit is to become 1 or 0).  Note that on return from the 
//    first call, CHANGE is set to 0, indicating that we begin with
//    the empty set.
//
{
  static int *a = NULL;
  int i;
  static int k = 0;
  static int n_save = -1;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "GRAY_NEXT - Fatal error!\n";
    cerr << "  Input value of N <= 0.\n";
    exit ( 1 );
  }

  if ( change == -n )
  {
    if ( a )
    {
      delete [] a;
    }

    a = new int[n];
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0;
    }

    n_save = n;
    k = 1;
    change = 0;

    return;
  }

  if ( n != n_save )
  {
    cerr << "\n";
    cerr << "GRAY_NEXT - Fatal error!\n";
    cerr << "  Input value of N has changed from definition value.\n";
    exit ( 1 );
  }
//
//  First determine WHICH item is to be changed.
//
  if ( ( k % 2 ) == 1 )
  {
    change = 1;
  }
  else
  {
    for ( i = 1; i <= n_save; i++ )
    {
      if ( a[i-1] == 1 )
      {
        change = i + 1;
        break;
      }
    }
  }
//
//  Take care of the terminal case CHANGE = N + 1.
//
  if ( change == n + 1 )
  {
    change = n;
  }
//
//  Now determine HOW the item is to be changed.
//
  if ( a[change-1] == 0 )
  {
    a[change-1] = 1;
  }
  else if ( a[change-1] == 1 )
  {
    a[change-1] = 0;
    change = -( change );
  }
//
//  Update the counter.
//
  k = k + 1;
//
//  If the output CHANGE = -N_SAVE, then we're done.
//
  if ( change == -n_save )
  {
    delete [] a;
    a = NULL;
    n_save = 0;
    k = 0;
  }

  return;
}
//****************************************************************************80

int gray_rank ( int gray )

//****************************************************************************80
//
//  Purpose:
//
//    GRAY_RANK ranks a Gray code.
//
//  Discussion:
//
//    Given the number GRAY, its ranking is the order in which it would be
//    visited in the Gray code ordering.  The Gray code ordering begins
//
//    Rank  Gray  Gray
//          (Dec) (Bin)
//
//       0     0  0000
//       1     1  0001
//       2     3  0011
//       3     2  0010
//       4     6  0110
//       5     7  0111
//       6     5  0101
//       7     4  0100
//       8    12  0110
//       etc
//
//   This routine is given a Gray code, and has to return the rank.
//
//  Example:
//
//    Gray  Gray  Rank
//    (Dec) (Bin)
//
//     0       0     0
//     1       1     1
//     2      10     3
//     3      11     2
//     4     100     7
//     5     101     6
//     6     110     4
//     7     111     5
//     8    1000    15
//     9    1001    14
//    10    1010    12
//    11    1011    13
//    12    1100     8
//    13    1101     9
//    14    1110    11
//    15    1111    10
//    16   10000    31
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int GRAY, the Gray code to be ranked.
//
//    Output, int GRAY_RANK, the rank of GRAY, and the integer whose Gray
//    code is GRAY.
//
{
  int i;
  int nbits = 32;
  int rank;

  rank = 0;

  if ( i4_btest ( gray, nbits - 1 ) )
  {
    rank = i4_bset ( rank, nbits - 1 );
  }

  for ( i = nbits-2; 0 <= i; i-- )
  {
    if ( i4_btest ( rank, i + 1 ) != i4_btest ( gray, i ) )
    {
      rank = i4_bset ( rank, i );
    }
  }
  return rank;
}
//****************************************************************************80

int gray_rank2 ( int gray )

//****************************************************************************80
//
//  Purpose:
//
//    GRAY_RANK2 ranks a Gray code.
//
//  Discussion:
//
//    In contrast to GRAY_RANK, this routine is entirely arithmetical,
//    and does not require access to bit testing and setting routines.
//
//
//    Given the number GRAY, its ranking is the order in which it would be
//    visited in the Gray code ordering.  The Gray code ordering begins
//
//    Rank  Gray  Gray
//          (Dec) (Bin)
//
//       0     0  0000
//       1     1  0001
//       2     3  0011
//       3     2  0010
//       4     6  0110
//       5     7  0111
//       6     5  0101
//       7     4  0100
//       8    12  0110
//       etc
//
//   This routine is given a Gray code, and has to return the rank.
//
//  Example:
//
//    Gray  Gray  Rank
//    (Dec) (Bin)
//
//     0       0     0
//     1       1     1
//     2      10     3
//     3      11     2
//     4     100     7
//     5     101     6
//     6     110     4
//     7     111     5
//     8    1000    15
//     9    1001    14
//    10    1010    12
//    11    1011    13
//    12    1100     8
//    13    1101     9
//    14    1110    11
//    15    1111    10
//    16   10000    31
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int GRAY, the Gray code to be ranked.
//
//    Output, int GRAY_RANK, the rank of GRAY, and the integer whose Gray
//    code is GRAY.
//
{
  int k;
  bool last;
  bool next;
  int rank;
  int two_k;

  if ( gray < 0 )
  {
    cerr << "\n";
    cerr << "GRAY_RANK2 - Fatal error!\n";
    cerr << "  Input value of GRAY < 0.\n";
    exit ( 1 );
  }

  if ( gray == 0 )
  {
    rank = 0;
    return rank;
  }
//
//  Find TWO_K, the largest power of 2 less than or equal to GRAY.
//
  k = 0;
  two_k = 1;
  while ( 2 * two_k <= gray )
  {
    two_k = two_k * 2;
    k = k + 1;
  }

  rank = two_k;
  last = true;
  gray = gray - two_k;

  while ( 0 < k )
  {
    two_k = two_k / 2;
    k = k - 1;

    next = ( two_k <= gray && gray < two_k * 2 );

    if ( next )
    {
      gray = gray - two_k;
    }

    if ( next != last )
    {
      rank = rank + two_k;
      last = true;
    }
    else
    {
      last = false;
    }
  }
  return rank;
}
//****************************************************************************80

int gray_unrank ( int rank )

//****************************************************************************80
//
//  Purpose:
//
//    GRAY_UNRANK unranks a Gray code.
//
//  Discussion:
//
//    The binary values of the Gray codes of successive integers differ in
//    just one bit.
//
//    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
//    Hamiltonian cycle on a graph of the cube in N dimensions.
//
//  Example:
//
//    Rank  Gray  Gray
//          (Dec) (Bin)
//
//     0     0       0
//     1     1       1
//     2     3      11
//     3     2      10
//     4     6     110
//     5     7     111
//     6     5     101
//     7     4     100
//     8    12    1100
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int RANK, the integer whose Gray code is desired.
//
//    Output, int GRAY_UNRANK, the Gray code of the given rank.
//
{
  int gray;
  int i;
  int nbits = 32;

  gray = 0;

  if ( i4_btest ( rank, nbits-1 ) )
  {
    gray = i4_bset ( gray, nbits-1 );
  }

  for ( i = nbits-2; 0 <= i; i-- )
  {
    if ( i4_btest ( rank, i+1 ) !=  i4_btest ( rank, i ) )
    {
      gray = i4_bset ( gray, i );
    }
  }
  return gray;
}
//****************************************************************************80

int gray_unrank2 ( int rank )

//****************************************************************************80
//
//  Purpose:
//
//    GRAY_UNRANK2 unranks a Gray code.
//
//  Discussion:
//
//    In contrast to GRAY_UNRANK, this routine is entirely arithmetical,
//    and does not require access to bit testing and setting routines.
//
//    The binary values of the Gray codes of successive integers differ in
//    just one bit.
//
//    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
//    Hamiltonian cycle on a graph of the cube in N dimensions.
//
//  Example:
//
//    Rank  Gray  Gray
//          (Dec) (Bin)
//
//     0     0       0
//     1     1       1
//     2     3      11
//     3     2      10
//     4     6     110
//     5     7     111
//     6     5     101
//     7     4     100
//     8    12    1100
//     9    14    1001
//    10    12    1010
//    11    13    1011
//    12     8    1100
//    13     9    1101
//    14    11    1110
//    15    10    1111
//    16    31   10000
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int RANK, the integer whose Gray code is desired.
//
//    Output, int GRAY_UNRANK2, the Gray code of the given rank.
//
{
  int gray;
  int k;
  bool last;
  bool next;
  int two_k;

  if ( rank <= 0 )
  {
    gray = 0;
    return gray;
  }

  k = 0;
  two_k = 1;
  while ( 2 * two_k <= rank )
  {
    two_k = two_k * 2;
    k = k + 1;
  }

  gray = two_k;
  rank = rank - two_k;
  next = true;

  while ( 0 < k )
  {
    two_k = two_k / 2;
    k = k - 1;

    last = next;
    next = ( two_k <= rank && rank <= two_k * 2 );

    if ( next != last )
    {
      gray = gray + two_k;
    }

    if ( next )
    {
      rank = rank - two_k;
    }
  }
  return gray;
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
      j = ( i4_huge ( ) + i4 ) + 1;
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
      value = ( i4_huge ( ) + i4 ) + 1;
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
      j = ( i4_huge ( ) + i4 ) + 1;
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
      value = - ( i4_huge ( ) - i4 ) - 1;
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
      j = ( i4_huge ( ) + i4 ) + 1;
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

int i4_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE computes the binomial coefficient C(N,K).
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in integer arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, are the values of N and K.
//
//    Output, int I4_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0;
  }
  else if ( mn == 0 )
  {
    value = 1;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( mx + i ) ) / i;
    }
  }
  return value;
}
//****************************************************************************80

void i4_factor ( int n, int maxfactor, int &nfactor, int factor[], 
  int exponent[], int &nleft )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTOR factors an I4 into prime factors.
//
//  Discussion:
//
//    N = NLEFT * Product ( 1 <= I <= NFACTOR ) FACTOR(I)^EXPONENT(I).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer to be factored.  N may be positive,
//    negative, or 0.
//
//    Input, int MAXFACTOR, the maximum number of prime factors for
//    which storage has been allocated.
//
//    Output, int &NFACTOR, the number of prime factors of N discovered
//    by the routine.
//
//    Output, int FACTOR[MAXFACTOR], the prime factors of N.
//
//    Output, int EXPONENT[MAXFACTOR].  EXPONENT(I) is the power of
//    the FACTOR(I) in the representation of N.
//
//    Output, int &NLEFT, the factor of N that the routine could not
//    divide out.  If NLEFT is 1, then N has been completely factored.
//    Otherwise, NLEFT represents factors of N involving large primes.
//
{
  int i;
  int maxprime;
  int p;

  nfactor = 0;

  for ( i = 0; i < maxfactor; i++ )
  {
    factor[i] = 0;
  }

  for ( i = 0; i < maxfactor; i++ )
  {
    exponent[i] = 0;
  }

  nleft = n;

  if ( n == 0 )
  {
    return;
  }

  if ( abs ( n ) == 1 )
  {
    nfactor = 1;
    factor[0] = 1;
    exponent[0] = 1;
    return;
  }
//
//  Find out how many primes we stored.
//
  maxprime = prime ( -1 );
//
//  Try dividing the remainder by each prime.
//
  for ( i = 1; i <= maxprime; i++ )
  {
    p = prime ( i );

    if ( abs ( nleft ) % p == 0 )
    {
      if ( nfactor < maxfactor )
      {
        nfactor = nfactor + 1;
        factor[nfactor-1] = p;
        exponent[nfactor-1] = 0;

        for ( ; ; )
        {
          exponent[nfactor-1] = exponent[nfactor-1] + 1;
          nleft = nleft / p;

          if ( abs ( nleft ) % p != 0 )
          {
            break;
          }
        }

        if ( abs ( nleft ) == 1 )
        {
          break;
        }
      }
    }
  }
  return;
}
//****************************************************************************80

int i4_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//    0 <= N <= 13 is required.
//
//    Output, int I4_FACTORIAL, the factorial of N.
//
{
  int i;
  int value;

  value = 1;

  if ( 13 < n ) 
  {
    value = - 1;
    cerr << "I4_FACTORIAL - Fatal error!\n";
    cerr << "  I4_FACTORIAL(N) cannot be computed as an integer\n";
    cerr << "  for 13 < N.\n";
    cerr << "  Input value N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 1; i <= n; i++ )
  {
    value = value * i;
  }

  return value;
}
//****************************************************************************80

int i4_gcd ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_GCD finds the greatest common divisor of I and J.
//
//  Discussion:
//
//    Only the absolute values of I and J are considered, so that the 
//    result is always nonnegative.
//
//    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
//
//    If I and J have no common factor, I4_GCD is returned as 1.
//
//    Otherwise, using the Euclidean algorithm, I4_GCD is the
//    largest common factor of I and J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, two numbers whose greatest common divisor
//    is desired.
//
//    Output, int I4_GCD, the greatest common divisor of I and J.
//
{
  int ip;
  int iq;
  int ir;
//
//  Return immediately if either I or J is zero.
//
  if ( i == 0 )
  {
    return i4_max ( 1, abs ( j ) );
  }
  else if ( j == 0 )
  {
    return i4_max ( 1, abs ( i ) );
  }
//
//  Set IP to the larger of I and J, IQ to the smaller.
//  This way, we can alter IP and IQ as we go.
//
  ip = i4_max ( abs ( i ), abs ( j ) );
  iq = i4_min ( abs ( i ), abs ( j ) );
//
//  Carry out the Euclidean algorithm.
//
  for ( ; ; )
  {
    ir = ip % iq;

    if ( ir == 0 )
    {
      break;
    }

    ip = iq;
    iq = ir;
  }
  return iq;
}
//****************************************************************************80

int i4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" integer value, usually the largest legal signed int.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int I4_HUGE, a "huge" integer.
//
{
  return 2147483647;
}
//****************************************************************************80

int i4_log_10 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_10 returns the whole part of the logarithm base 10 of an integer.
//
//  Discussion:
//
//    It should be the case that 10^I_LOG_10(I) <= |I| < 10^(I_LOG_10(I)+1).
//    (except for I = 0).
//
//    The number of decimal digits in I is I4_LOG_10(I) + 1.
//
//  Example:
//
//        I    I4_LOG_10(I)
//
//        0     0
//        1     0
//        2     0
//
//        9     0
//       10     1
//       11     1
//
//       99     1
//      100     2
//      101     2
//
//      999     2
//     1000     3
//     1001     3
//
//     9999     3
//    10000     4
//    10001     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer.
//
//    Output, int I4_LOG_10, the whole part of the logarithm of abs ( I ).
//
{
  int ten_pow;
  int value;

  i = abs ( i );

  ten_pow = 10;
  value = 0;

  while ( ten_pow <= i )
  {
    ten_pow = ten_pow * 10;
    value = value + 1;
  }

  return value;
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
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of i1 and i2.
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
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
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

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If 
//      NREM = I4_MODP ( I, J ) 
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Example:
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
// 
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is 
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_moebius ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MOEBIUS returns the value of MU(N), the Moebius function of N.
//
//  Discussion:
//
//    MU(N) is defined as follows:
//
//      MU(N) = 1 if N = 1;
//              0 if N is divisible by the square of a prime;
//              (-1)**K, if N is the product of K distinct primes.
//
//    The Moebius function MU(D) is related to Euler's totient 
//    function PHI(N):
//
//      PHI(N) = sum ( D divides N ) MU(D) * ( N / D ).
//
//    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
//    if N is a square, cube, etc.
//
//  Example:
//
//     N  MU(N)
//
//     1    1
//     2   -1
//     3   -1
//     4    0
//     5   -1
//     6    1
//     7   -1
//     8    0
//     9    0
//    10    1
//    11   -1
//    12    0
//    13   -1
//    14    1
//    15    1
//    16    0
//    17   -1
//    18    0
//    19   -1
//    20    0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the value to be analyzed.
//
//    Output, int I4_MOEBIUS, the value of MU(N).
//    If N is less than or equal to 0, MU will be returned as -2.
//    If there was not enough internal space for factoring, MU
//    is returned as -3.
//
{
# define FACTOR_MAX 20

  int exponent[FACTOR_MAX];
  int factor[FACTOR_MAX];
  int i;
  int mu;
  int nfactor;
  int nleft;

  if ( n <= 0 )
  {
    mu = -2;
    return mu;
  }

  if ( n == 1 )
  {
    mu = 1;
    return mu;
  }
//
//  Factor N.
//
  i4_factor ( n, FACTOR_MAX, nfactor, factor, exponent, nleft );

  if ( nleft != 1 )
  {
    cerr << "\n";
    cerr << "I4_MOEBIUS - Fatal error!\n";
    cerr << "  Not enough factorization space.\n";
    mu = -3;
    exit ( 1 );
  }

  mu = 1;

  for ( i = 0; i < nfactor; i++ )
  {
    mu = -mu;

    if ( 1 < exponent[i] )
    {
      mu = 0;
      return mu;
    }
  }

  return mu;
# undef FACTOR_MAX
}
//****************************************************************************80

void i4_partition_conj ( int n, int iarray1[], int mult1[], int npart1, 
  int iarray2[], int mult2[], int &npart2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITION_CONJ computes the conjugate of a partition.
//
//  Discussion:
//
//    A partition of an integer N is a set of positive integers which
//    add up to N.  The conjugate of a partition P1 of N is another partition
//    P2 of N obtained in the following way:
//
//      The first element of P2 is the number of parts of P1 greater than
//      or equal to 1.
//
//      The K-th element of P2 is the number of parts of P1 greater than
//      or equal to K.
//
//    Clearly, P2 will have no more than N elements; it may be surprising
//    to find that P2 is guaranteed to be a partition of N.  However, if
//    we symbolize the initial partition P1 by rows of X's, then we can
//    see that P2 is simply produced by grouping by columns:
//
//        6 3 2 2 1
//      5 X X X X X
//      4 X X X X
//      2 X X
//      1 X
//      1 X
//      1 X
//
//  Example:
//
//    14 = 5 + 4 + 2 + 1 + 1 + 1
//
//    The conjugate partition is:
//
//    14 = 6 + 3 + 2 + 2 + 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int N, the integer to be partitioned.
//
//    Input, int IARRAY1[NPART1], contains the parts of
//    the partition.  The value of N is represented by
//
//      sum ( 1 <= I <= NPART1 ) MULT1(I) * IARRAY1(I).
//
//    Input, int MULT1[NPART1], counts the multiplicity of
//    the parts of the partition.  MULT1(I) is the multiplicity
//    of the part IARRAY1(I), for 1 <= I <= NPART1.
//
//    Input, int NPART1, the number of "parts" in the partition.
//
//    Output, int IARRAY2[N], contains the parts of
//    the conjugate partition in entries 1 through NPART2.
//
//    Output, int MULT2[N], counts the multiplicity of
//    the parts of the conjugate partition in entries 1 through NPART2.
//
//    Output, int &NPART2, the number of "parts" in the conjugate partition.
//
{
  int i;
  int itemp;
  int itest;

  for ( i = 0; i < n; i++ )
  {
    iarray2[i] = 0;
  }
  for ( i = 0; i < n; i++ )
  {
    mult2[i] = 0;
  }
  npart2 = 0;

  itest = 0;

  for ( ; ; )
  {
    itest = itest + 1;

    itemp = 0;

    for ( i = 0; i < npart1; i++ )
    {
      if ( itest <= iarray1[i] )
      {
        itemp = itemp + mult1[i];
      }
    }

    if ( itemp <= 0 )
    {
      break;
    }

    if ( 0 < npart2 )
    {
      if ( itemp == iarray2[npart2-1] )
      {
        mult2[npart2-1] = mult2[npart2-1] + 1;
      }
      else
      {
        npart2 = npart2 + 1;
        iarray2[npart2-1] = itemp;
        mult2[npart2-1] = 1;
      }
    }
    else
    {
      npart2 = npart2 + 1;
      iarray2[npart2-1] = itemp;
      mult2[npart2-1] = 1;
    }
  }
  return;
}
//****************************************************************************80

void i4_partition_count ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITION_COUNT computes the number of partitions of an integer.
//
//  Discussion:
//
//    Partition numbers are difficult to compute.  This routine uses
//    Euler's method, which observes that:
//
//      P(0) = 1
//      P(N) =   P(N-1)  + P(N-2)
//             - P(N-5)  - P(N-7)
//             + P(N-12) + P(N-15)
//             - ...
//
//      where the numbers 1, 2, 5, 7, ... to be subtracted from N in the
//      indices are the successive pentagonal numbers, (with both positive 
//      and negative indices) with the summation stopping when a negative 
//      index is reached.
//
//  First values:
//
//    N   P
//
//    0   1
//    1   1
//    2   2
//    3   3
//    4   5
//    5   7
//    6  11
//    7  15
//    8  22
//    9  30
//   10  42
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    John Conway, Richard Guy,
//    The Book of Numbers,
//    Springer Verlag, 1996, page 95.
//
//  Parameters:
//
//    Input, int N, the index of the highest partition number desired.
//
//    Output, int P[N+1], the partition numbers.
//
{
  int i;
  int j;
  int pj;
  int sgn;

  p[0] = 1;

  for ( i = 1; i <= n; i++ )
  {
    p[i] = 0;

    j = 0;
    sgn = 1;

    for ( ; ; )
    {
      j = j + 1;
      pj = pent_enum ( j );

      if ( i < pj )
      {
        break;
      }

      p[i] = p[i] + sgn * p[i-pj];
      sgn = -sgn;
    }

    j = 0;
    sgn = 1;

    for ( ; ; )
    {
      j = j - 1;
      pj = pent_enum ( j );

      if ( i < pj )
      {
        break;
      }

      p[i] = p[i] + sgn * p[i-pj];
      sgn = -sgn;

    }
  }

  return;
}
//****************************************************************************80

int *i4_partition_count2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITION_COUNT2 computes the number of partitions of an integer.
//
//  First values:
//
//    N   P
//
//    0   1
//    1   1
//    2   2
//    3   3
//    4   5
//    5   7
//    6  11
//    7  15
//    8  22
//    9  30
//   10  42
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the largest integer to be considered.
//
//    Output, int I4_PARTITION_COUNT2[0:N], the partition numbers.
//
{
  int i;
  int j;
  int *p;
  int s;
  int t;
  int total;

  if ( n < 0 )
  {
    return NULL;
  }

  p = new int[n+1];

  p[0] = 1;

  if ( n < 1 )
  {
    return p;
  }

  p[1] = 1;

  for ( i = 2; i <= n; i++ )
  {
    total = 0;

    for ( t = 1; t <= i; t++ )
    {
      s = 0;
      j = i;

      for ( ; ; )
      {
        j = j - t;

        if ( 0 < j )
        {
          s = s + p[j];
        }
        else
        {
          if ( j == 0 )
          {
            s = s + 1;
          }
          break;
        }
      }
      total = total + s * t;
    }
    p[i] = ( int ) ( total / i );
  }
  return p;
}
//****************************************************************************80

void i4_partition_count_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITION_COUNT_VALUES returns some values of the int partition count.
//
//  Discussion:
//
//    A partition of an integer I is a representation of the integer
//    as the sum of nonzero positive integers.  The order of the summands
//    does not matter.  Thus, the number 5 has the following partitions
//    and no more:
//
//    5 = 5
//      = 4 + 1 
//      = 3 + 2 
//      = 3 + 1 + 1 
//      = 2 + 2 + 1 
//      = 2 + 1 + 1 + 1 
//      = 1 + 1 + 1 + 1 + 1
//
//    so the number of partitions of 5 is 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int &N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When 
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int &N, the integer.
//
//    Output, int &C, the number of partitions of the integer.
//
{
# define N_MAX 21

  int c_vec[N_MAX] = {
      1,
      1,   2,   3,   5,   7,  11,  15,  22,  30,  42, 
     56,  77, 101, 135, 176, 231, 297, 385, 490, 627 };
  int n_vec[N_MAX] = {
     0,  
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  if ( N_MAX <= n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data];
    c = c_vec[n_data];
    n_data = n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void i4_partition_next ( bool &done, int a[], int mult[], int n, int &npart )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITION_NEXT generates the partitions of an integer, one at a time.
//
//  Discussion:
//
//    The number of partitions of N is:
//
//      1     1
//      2     2
//      3     3
//      4     5
//      5     7
//      6    11
//      7    15
//      8    22
//      9    30
//     10    42
//     11    56
//     12    77
//     13   101
//     14   135
//     15   176
//     16   231
//     17   297
//     18   385
//     19   490
//     20   627
//     21   792
//     22  1002
//     23  1255
//     24  1575
//     25  1958
//     26  2436
//     27  3010
//     28  3718
//     29  4565
//     30  5604
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input/output, bool &DONE.
//    On first call, the user should set DONE to TRUE to signal
//    that the program should initialize data.
//    On each return, the programs sets DONE to FALSE if it
//    has another partition to return.  If the program returns
//    with DONE TRUE, then there are no more partitions.
//
//    Output, int A[N].  A contains the parts of
//    the partition.  The value of N is represented by
//      N = sum ( 1 <= I <= NPART ) MULT(I) * A(I).
//
//    Output, int MULT[N].  MULT counts the multiplicity of
//    the parts of the partition.  MULT(I) is the multiplicity
//    of the part A(I), for 1 <= I <= NPART.
//
//    Input, int N, the integer to be partitioned.
//
//    Output, int &NPART, the number of "parts" in the partition.
//
{
  int i;
  int is;
  int iu;
  int iv;
  int iw;
  int k;
  int k1;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "I4_PARTITION_NEXT - Fatal error!\n";
    cerr << "  N must be positive.\n";
    cerr << "  The input value of N was " << n << "\n";
    exit ( 1 );
  }

  if ( done )
  {
    a[0] = n;
    for ( i = 1; i < n; i++ )
    {
      a[i] = 0;
    }

    mult[0] = 1;
    for ( i = 1; i < n; i++ )
    {
      mult[i] = 0;
    }
    npart = 1;
    done = false;
  }
  else
  {
    if ( 1 < a[(npart)-1] || 1 < npart )
    {
      done = false;

      if ( a[(npart)-1] == 1 )
      {
        is = a[npart-2] + mult[npart-1];
        k = npart - 1;
      }
      else
      {
        is = a[npart-1];
        k = npart;
      }

      iw = a[k-1] - 1;
      iu = is / iw;
      iv = is % iw;
      mult[k-1] = mult[k-1] - 1;

      if ( mult[k-1] == 0 )
      {
        k1 = k;
      }
      else
      {
        k1 = k + 1;
      }

      mult[k1-1] = iu;
      a[k1-1] = iw;

      if ( iv == 0 )
      {
        npart = k1;
      }
      else
      {
        mult[k1] = 1;
        a[k1] = iv;
        npart = k1 + 1;
      }
    }
    else
    {
      done = true;
    }
  }
  return;
}
//****************************************************************************80

void i4_partition_next2 ( int n, int a[], int mult[], int &npart, bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITION_NEXT2 computes the partitions of an integer one at a time.
//
//  Discussion:
//
//    Unlike compositions, order is not important in a partition.  Thus
//    the sequences 3+2+1 and 1+2+3 represent distinct compositions, but
//    not distinct partitions.  Also 0 is never returned as one of the
//    elements of the partition.
//
//  Example:
//
//    Sample partitions of 6 include:
//
//      6 = 4+1+1 = 3+2+1 = 2+2+2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the integer whose partitions are desired.
//
//    Output, int A[N].  A(I) is the I-th distinct part
//    of the partition, for I = 1, NPART.  Note that if a certain number
//    shows up several times in the partition, it is listed only
//    once in A, and its multiplicity is counted in MULT.
//
//    Output, int MULT[N].  MULT(I) is the multiplicity of A(I)
//    in the partition, for I = 1, NPART; that is, the number of repeated
//    times that A(I) is used in the partition.
//
//    Output, int &NPART, the number of distinct, nonzero parts in the
//    output partition.
//
//    Input/output, bool &MORE.  Set MORE = FALSE on first call.  It
//    will be reset TRUE on return with the first partition.
//    Keep calling for more partitions until MORE
//    is returned FALSE
//
{
  int iff;
  int is;
  int isum;
  static int nlast = 0;
//
//  On the first call, set NLAST to 0.
//
  if ( !more )
  {
    nlast = 0;
  }

  if ( n != nlast || ( !more ) )
  {
    nlast = n;
    npart = 1;
    a[npart-1] = n;
    mult[npart-1] = 1;
    more = mult[npart-1] != n;
    return;
  }

  isum = 1;

  if ( a[npart-1] <= 1 )
  {
    isum = mult[npart-1] + 1;
    npart = npart - 1;
  }

  iff = a[npart-1] - 1;

  if ( mult[npart-1] != 1 )
  {
    mult[npart-1] = mult[npart-1] - 1;
    npart = npart + 1;
  }

  a[npart-1] = iff;
  mult[npart-1] = 1 + ( isum / iff );
  is = isum % iff;

  if ( 0 < is )
  {
    npart = npart + 1;
    a[npart-1] = is;
    mult[npart-1] = 1;
  }
//
//  There are more partitions, as long as we haven't just computed
//  the last one, which is N copies of 1.
//
  more = mult[npart-1] != n;

  return;
}
//****************************************************************************80

void i4_partition_print ( int n, int npart, int a[], int mult[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITION_PRINT prints a partition of an integer.
//
//  Discussion:
//
//    A partition of an int N is a representation of the integer as
//    the sum of nonzero integers:
//
//      N = A1 + A2 + A3 + ...
//
//    It is standard practice to gather together all the values that 
//    are equal, and replace them in the sum by a single term, multiplied
//    by its "multiplicity":
//
//      N = M1 * A1 + M2 * A2 + ... + M(NPART) * A(NPART)
//    
//    In this representation, every A is a unique positive number, and 
//    no M is zero (or negative).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer to be partitioned.
//
//    Input, int NPART, the number of "parts" in the partition.
//
//    Input, int A[NPART], the parts of the partition.  
//
//    Input, int MULT[NPART], the multiplicities of the parts.
//
{
  int i;

  cout << "  " << n << " = ";

  for ( i = 0; i < npart; i++ )
  {
    if ( 0 < i )
    {
      cout << " + ";
    }
    cout << mult[i] << " * " << a[i];
  }
  cout << "\n";

  return;
}
//****************************************************************************80

void i4_partition_random ( int n, int table[], int &seed, int a[], int mult[], 
  int &npart )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITION_RANDOM selects a random partition of the int N.
//
//  Discussion:
//
//    Note that some elements of the partition may be 0.  The partition is
//    returned as (MULT(I),I), with NPART nonzero entries in MULT, and
//
//      N = sum ( 1 <= I <= N ) MULT(I) * I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer to be partitioned.
//
//    Input, int TABLE[N], the number of partitions of each integer 
//    from 1 to N.  This table may be computed by I4_PARTITION_COUNT2.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int A[N], contains in A(1:NPART) the parts of the partition.
//
//    Output, int MULT[N], contains in MULT(1:NPART) the multiplicity
//    of the parts.
//
//    Output, int &NPART, the number of parts in the partition chosen,
//    that is, the number of integers I with nonzero multiplicity MULT(I).
//
{
  int i;
  int i1;
  int id;
  int j;
  int m;
  double z;

  m = n;
  npart = 0;
  for ( i = 0; i < n; i++ )
  {
    mult[i] = 0;
  }

  while ( 0 < m )
  {
    z = r8_uniform_01 ( seed );
    z = m * table[m-1] * z;
    id = 1;
    i1 = m;
    j = 0;

    for ( ; ; )
    {
      j = j + 1;
      i1 = i1 - id;

      if ( i1 < 0 )
      {
        id = id + 1;
        i1 = m;
        j = 0;
        continue;
      }

      if ( i1 == 0 )
      {
        z = z - id;
        if ( 0.0 < z )
        {
          id = id + 1;
          i1 = m;
          j = 0;
          continue;
        }
        else
        {
          break;
        }
      }

      if ( 0 < i1 )
      {
        z = z - id * table[i1-1];
        if ( z <= 0.0 )
        {
          break;
        }
      }
    }

    mult[id-1] = mult[id-1] + j;
    npart = npart + j;
    m = i1;
  }
//
//  Reformulate the partition in the standard form.
//  NPART is the number of distinct parts.
//
  npart = 0;

  for ( i = 1; i <= n; i++ )
  {
    if ( mult[i-1] != 0 )
    {
      npart = npart + 1;
      a[npart-1] = i;
      mult[npart-1] = mult[i-1];
    }
  }

  for ( i = npart + 1; i <= n; i++ )
  {
    mult[i-1] = 0;
  }

  return;
}
//****************************************************************************80

void i4_partitions_next ( int s, int m[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITIONS_NEXT: next partition into S parts.
//
//  Discussion:
//
//    This function generates, one at a time, entries from the list of
//    nondecreasing partitions of the integers into S or fewer parts.
//
//    The list is ordered first by the integer that is partitioned
//    (the sum of the entries), and second by decreasing lexical order
//    in the partition vectors.
//
//    The first value returned is the only such partition of 0.
//
//    Next comes the only partition of 1.
//
//    There follow two partitions of 2, and so on.
//
//    Typical use of this function begins with an initialization call,
//    and then repeated calls in which the output from the previous call
//    is used as input to the next call:
//
//    m = [ 0, 0, 0 ];
//
//    while ( condition )
//      m = i4_partitions_next ( s, m );
//    end
//
//  Example:
//
//    S = 3
//
//    P  D    M
//    _  _  _____
//    1  0  0 0 0
//    2  1  1 0 0
//    3  2  2 0 0
//    4  2  1 1 0
//    5  3  3 0 0
//    6  3  2 1 0
//    7  3  1 1 1
//    8  4  4 0 0
//    9  4  3 1 0
//   10  4  2 2 0
//   11  4  2 1 1
//   12  5  5 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2010
//
//  Author:
//
//    Original MATLAB version by Alan Genz.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int S, the number of items in the partition.
//
//    Input/output, int M[S].  On input, the current partition.  
//    On first call, this should be a nondecreasing partition.  Thereafter, it 
//    should be the output partition from the previous call.  On output, the
//    next partition.
//
{
  int i;
  int j;
  int msum;

  msum = m[0];

  for ( i = 1; i < s; i++ )
  {
    msum = msum + m[i];

    if ( m[0] <= m[i] + 1 )
    {
      m[i] = 0;
    }
    else
    {
      m[0] = msum - i * ( m[i] + 1 );
      for ( j = 1; j <= i; j++ )
      {
        m[j] = m[i] + 1;
      }
      return;
    }
  }
//
//  If we failed to find a suitable index I, put
//  the entire sum into M(1), increment by 1, and
//  prepare to partition the next integer.
//
  m[0] = msum + 1;

  return;
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
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
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
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
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

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Discussion:
//
//    The sign of 0 and all positive integers is taken to be +1.
//    The sign of all negative integers is -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
  if ( i < 0 ) 
  {
    return (-1);
  }
  else
  {
    return 1;
  }

}
//****************************************************************************80

void i4_sqrt ( int n, int &q, int &r )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SQRT finds the integer square root of N by solving N = Q^2 + R.
//
//  Discussion:
//
//    The integer square root of N is an integer Q such that
//    Q**2 <= N but N < (Q+1)^2.
//
//    A simpler calculation would be something like
//
//      Q = INT ( SQRT ( REAL ( N ) ) )
//
//    but this calculation has the virtue of using only integer arithmetic.
//
//    To avoid the tedium of worrying about negative arguments, the routine
//    automatically considers the absolute value of the argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2003
//
//  Author:
//
//   John Burkardt
//
//  Reference:
//
//    Mark Herkommer,
//    Number Theory, A Programmer's Guide,
//    McGraw Hill, 1999, pages 294-307.
//
//  Parameters:
//
//    Input, int N, the number whose integer square root is desired.
//    Actually, only the absolute value of N is considered.
//
//    Output, int &Q, &R, the integer square root, and positive remainder,
//    of N.
//
{
  n = abs ( n );

  q = n;

  if ( 0 < n )
  {
    while ( ( n / q ) < q )
    {
      q = ( q + ( n / q ) ) / 2;
    }
  }

  r = n - q * q;

  return;
}
//****************************************************************************80

void i4_sqrt_cf ( int n, int max_term, int &n_term, int b[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SQRT_CF finds the continued fraction representation of a square root of an integer.
//
//  Discussion:
//
//    The continued fraction representation of the square root of an integer
//    has the form
//
//      [ B0, (B1, B2, B3, ..., BM), ... ]
//
//    where
//
//      B0 = int ( sqrt ( real ( N ) ) )
//      BM = 2 * B0
//      the sequence ( B1, B2, B3, ..., BM ) repeats in the representation.
//      the value M is termed the period of the representation.
//
//  Example:
//
//     N  Period  Continued Fraction
//
//     2       1  [ 1, 2, 2, 2, ... ]
//     3       2  [ 1, 1, 2, 1, 2, 1, 2... ]
//     4       0  [ 2 ]
//     5       1  [ 2, 4, 4, 4, ... ]
//     6       2  [ 2, 2, 4, 2, 4, 2, 4, ... ]
//     7       4  [ 2, 1, 1, 1, 4, 1, 1, 4, 1, 1, 4... ]
//     8       2  [ 2, 1, 4, 1, 4, 1, 4, 1, 4, ... ]
//     9       0  [ 3 ]
//    10       1  [ 3, 6, 6, 6, ... ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//   John Burkardt
//
//  Reference:
//
//    Mark Herkommer,
//    Number Theory, A Programmer's Guide,
//    McGraw Hill, 1999, pages 294-307.
//
//  Parameters:
//
//    Input, int N, the number whose continued fraction square root
//    is desired.
//
//    Input, int MAX_TERM, the maximum number of terms that may
//    be computed.
//
//    Output, int &N_TERM, the number of terms computed beyond the
//    0 term.  The routine should stop if it detects that the period
//    has been reached.
//
//    Output, int B[MAX_TERM+1], contains the continued fraction
//    coefficients for indices 0 through N_TERM.
//
{
  int k;
  int p;
  int q;
  int r;
  int s;

  n_term = 0;

  i4_sqrt ( n, s, r );
  b[0] = s;

  if ( 0 < r )
  {
    p = 0;
    q = 1;

    for ( ; ; )
    {
      p = b[n_term] * q - p;
      q = ( n - p * p ) / q;

      if ( max_term <= n_term )
      {
        return;
      }

      n_term = n_term + 1;
      b[n_term] = ( p + s ) / q;

      if ( q == 1 )
      {
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

void i4_swap ( int &i, int &j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &I, &J.  On output, the values of I and
//    J have been interchanged.
//
{
  int k;

  k = i;
  i = j;
  j = k;
 
  return;
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

void i4_to_chinese ( int j, int n, int m[], int r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_CHINESE converts an I4 to its Chinese remainder form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int J, the integer to be converted.
//
//    Input, int N, the number of moduluses.
//
//    Input, int M[N], the moduluses.  These should be positive
//    and pairwise prime.
//
//    Output, int R[N], the Chinese remainder representation of the integer.
//
{
  bool error;
  int i;

  error = chinese_check ( n, m );

  if ( error )
  {
    cerr << "\n";
    cerr << "I4_TO_CHINESE - Fatal error!\n";
    cerr << "  The moduluses are not legal.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    r[i] = i4_modp ( j, m[i] );
  }

  return;
}
//****************************************************************************80

void i4_to_dvec ( int i4, int n, int dvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_DVEC makes a signed decimal vector from an I4.
//
//  Discussion:
//
//    A DVEC is an integer vector of decimal digits, intended to
//    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
//    is the coefficient of 10**(N-2), and DVEC(N) contains sign
//    information.  It is 0 if the number is positive, and 9 if
//    the number is negative.
//
//    Negative values have a ten's complement operation applied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2006
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
//    Output, int DVEC[N], the signed decimal representation.
//
{
  int base = 10;
  int i;
  int i4_copy;

  i4_copy = abs ( i4 );

  for ( i = 0; i < n-1; i++ )
  {
    dvec[i] = i4_copy % base;

    i4_copy = i4_copy / base;
  }

  dvec[n-1] = 0;

  if ( i4 < 0 )
  {
    dvec_complementx ( n, dvec, dvec );
  }

  return;
}
//****************************************************************************80

void i4_to_i4poly ( int intval, int base, int degree_max, int &degree, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_I4POLY converts an I4 to an integer polynomial in a given base.
//
//  Example:
//
//    INTVAL  BASE  Degree     A (in reverse order!)
//
//         1     2       0     1
//         6     2       2     1  1  0
//        23     2       5     1  0  1  1  1
//        23     3       3     2  1  2
//        23     4       3     1  1  3
//        23     5       2     4  3
//        23     6       2     3  5
//        23    23       1     1  0
//        23    24       0    23
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int INTVAL, an integer to be converted.
//
//    Input, int BASE, the base, which should be greater than 1.
//
//    Input, int DEGREE_MAX, the maximum degree.
//
//    Output, int &DEGREE, the degree of the polynomial.
//
//    Output, int A[DEGREE_MAX+1], contains the coefficients
//    of the polynomial expansion of INTVAL in base BASE.
//
{
  int i;
  int j;

  for ( i = 0; i <= degree_max; i++ )
  {
    a[i] = 0;
  }

  j = abs ( intval );

  degree = 0;

  a[degree] = j % base;

  j = j - a[degree];
  j = j / base;

  while ( 0 < j )
  {
    degree = degree + 1;

    if ( degree < degree_max )
    {
      a[degree] = j % base;
    }

    j = j - a[degree];
    j = j / base;
  }

  if ( intval < 0 )
  {
    for ( i = 0; i <= degree_max; i++ )
    {
      a[i] = -a[i];
    }
  }

  return;
}
//****************************************************************************80

double i4_to_van_der_corput ( int seed, int base )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_VAN_DER_CORPUT computes an element of a van der Corput sequence.
//
//  Discussion:
//
//    The van der Corput sequence is often used to generate a "subrandom"
//    sequence of points which have a better covering property
//    than pseudorandom points.
//
//    The van der Corput sequence generates a sequence of points in [0,1]
//    which (theoretically) never repeats.  Except for SEED = 0, the
//    elements of the van der Corput sequence are strictly between 0 and 1.
//
//    The van der Corput sequence writes an integer in a given base B,
//    and then its digits are "reflected" about the decimal point.
//    This maps the numbers from 1 to N into a set of numbers in [0,1],
//    which are especially nicely distributed if N is one less
//    than a power of the base.
//
//    Hammersley suggested generating a set of N nicely distributed
//    points in two dimensions by setting the first component of the
//    Ith point to I/N, and the second to the van der Corput 
//    value of I in base 2.  
//
//    Halton suggested that in many cases, you might not know the number 
//    of points you were generating, so Hammersley's formulation was
//    not ideal.  Instead, he suggested that to generated a nicely
//    distributed sequence of points in M dimensions, you simply
//    choose the first M primes, P(1:M), and then for the J-th component of
//    the I-th point in the sequence, you compute the van der Corput
//    value of I in base P(J).
//
//    Thus, to generate a Halton sequence in a 2 dimensional space,
//    it is typical practice to generate a pair of van der Corput sequences,
//    the first with prime base 2, the second with prime base 3.
//    Similarly, by using the first K primes, a suitable sequence
//    in K-dimensional space can be generated.
//
//    The generation is quite simple.  Given an integer SEED, the expansion
//    of SEED in base BASE is generated.  Then, essentially, the result R
//    is generated by writing a decimal point followed by the digits of
//    the expansion of SEED, in reverse order.  This decimal value is actually
//    still in base BASE, so it must be properly interpreted to generate
//    a usable value.
//
//  Example:
//
//    BASE = 2
//
//    SEED     SEED      van der Corput
//    decimal  binary    binary   decimal
//    -------  ------    ------   -------
//        0  =     0  =>  .0     = 0.0
//        1  =     1  =>  .1     = 0.5
//        2  =    10  =>  .01    = 0.25
//        3  =    11  =>  .11    = 0.75
//        4  =   100  =>  .001   = 0.125
//        5  =   101  =>  .101   = 0.625
//        6  =   110  =>  .011   = 0.375
//        7  =   111  =>  .111   = 0.875
//        8  =  1000  =>  .0001  = 0.0625
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    John Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, 1960, pages 84-90.
// 
//    John Hammersley,
//    Monte Carlo methods for solving multivariable problems,
//    Proceedings of the New York Academy of Science,
//    Volume 86, 1960, pages 844-874.
//
//    Johannes van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Input, int SEED, the index of the desired element.
//    SEED should be nonnegative.
//    SEED = 0 is allowed, and returns R = 0.
//
//    Input, int BASE, the van der Corput base, which is usually 
//    a prime number.  BASE must be greater than 1.
//
//    Output, double VAN_DER_CORPUT, the SEED-th element of the van 
//    der Corput sequence for base BASE.
//
{
  double base_inv;
  int digit;
  double r;

  if ( base <= 1 )
  {
    cerr << "\n";
    cerr << "I4_TO_VAN_DER_CORPUT - Fatal error!\n";
    cerr << "  The input base BASE is <= 1!\n";
    cerr << "  BASE = " << base << "\n";
    exit ( 1 );
  }

  if ( seed < 0 ) 
  {
    cerr << "\n";
    cerr << "I4_TO_VAN_DER_CORPUT - Fatal error!\n";
    cerr << "  SEED < 0.";
    cerr << "  SEED = " << seed << "\n";
    exit ( 1 );
  }

  r = 0.0;

  base_inv = 1.0 / ( ( double ) base );

  while ( seed != 0 )
  {
    digit = seed % base;
    r = r + ( ( double ) digit ) * base_inv;
    base_inv = base_inv / ( ( double ) base );
    seed = seed / base;
  }

  return r;
}
//****************************************************************************80

int i4_uniform ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
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
//    12 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
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
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 ) 
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80

void i4mat_01_rowcolsum ( int m, int n, int r[], int c[], int a[], bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_01_ROWCOLSUM creates a 0/1 I4MAT with given row and column sums.
//
//  Discussion:
//
//    Given an M vector R and N vector C, there may exist one or more
//    M by N matrices with entries that are 0 or 1, whose row sums are R
//    and column sums are C.
//
//    For convenience, this routine requires that the entries of R and C
//    be given in nonincreasing order.
//
//    There are several requirements on R and C.  The simple requirements
//    are that the entries of R and C must be nonnegative, that the entries
//    of R must each be no greater than N, and those of C no greater than M,
//    and that the sum of the entries of R must equal the sum of the entries 
//    of C.
//
//    The final technical requirement is that if we form R*, the conjugate
//    partition of R, then C is majorized by R*, that is, that every partial
//    sum from 1 to K of the entries of C is no bigger than the sum of the same
//    entries of R*, for every K from 1 to N.
//
//    Given these conditions on R and C, there is at least one 0/1 matrix
//    with the given row and column sums.
//
//    The conjugate partition of R is constructed as follows:
//      R*(1) is the number of entries of R that are 1 or greater.
//      R*(2) is the number of entries of R that are 2 or greater.
//      ...
//      R*(N) is the number of entries of R that are N (can't be greater).
//
//  Example:
//
//    M = N = 5
//    R = ( 3, 2, 2, 1, 1 )
//    C = ( 2, 2, 2, 2, 1 )
//
//    A =
//      1 0 1 0 1
//      1 0 0 1 0
//      0 1 0 1 0
//      0 1 0 0 0
//      0 0 1 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack van Lint, Richard Wilson,
//    A Course in Combinatorics,
//    Oxford, 1992, pages 148-156.
//
//    James Sandeson,
//    Testing Ecobool Patterns,
//    American Scientist,
//    Volume 88, July-August 2000, pages 332-339.
//
//    Ian Saunders,
//    Algorithm AS 205,
//    Enumeration of R x C Tables with Repeated Row Totals,
//    Applied Statistics,
//    Volume 33, Number 3, pages 340-352, 1984.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, int R[M], C[N], the row and column sums desired for the array.
//    Both vectors must be arranged in descending order.
//    The elements of R must be between 0 and N.
//    The elements of C must be between 0 and M.
//
//    Output, int A[M*N], the M by N matrix with the given row and
//    column sums.
//    Each entry of A is 0 or 1.
//
//    Output, bool &ERROR, is true if an error occurred.
//
{
  int c_sum;
  int i;
  int j;
  int k;
  int *r_conj;
  int r_sum;
  int *r2;
//
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i+j*m] = 0;
    }
  }
//
//  Check conditions.
//
  error = false;

  if ( i4vec_sum ( m, r ) != i4vec_sum ( n, c ) )
  {
    cerr << "\n";
    cerr << "I4MAT_01_ROWCOLSUM - Fatal error!\n";
    cerr << "  Row sums R and column sums C don't have the same sum!\n";
    error = true;
    exit ( 1 );
  }

  if ( !i4vec_descends ( m, r ) )
  {
    cerr << "\n";
    cerr << "I4MAT_01_ROWCOLSUM - Fatal error!\n";
    cerr << "  Row sum vector R is not descending!\n";
    error = true;
    exit ( 1 );
  }

  if ( n < r[0] || r[m-1] < 0 )
  {
    error = true;
    return;
  }

  if ( !i4vec_descends ( n, c ) )
  {
    cerr << "\n";
    cerr << "I4MAT_01_ROWCOLSUM - Fatal error!\n";
    cerr << "  Column sum vector C is not descending!\n";
    error = true;
    exit ( 1 );
  }

  if ( m < c[0] || c[n-1] < 0 )
  {
    error = true;
    return;
  }
//
//  Compute the conjugate of R.
//
  r_conj = new int[n];

  for ( i = 0; i < n; i++ )
  {
    r_conj[i] = 0;
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < r[i]; j++ )
    {
      r_conj[j] = r_conj[j] + 1;
    }
  }
//
//  C must be majorized by R_CONJ.
//
  r_sum = 0;
  c_sum = 0;
  for ( i = 0; i < n; i++ )
  {
    r_sum = r_sum + r_conj[i];
    c_sum = c_sum + c[i];
    if ( r_sum < c_sum )
    {
      error = true;
      return;
    }
  }
  delete [] r_conj;

  if ( error )
  {
    return;
  }
  r2 = new int[m];

//
//  We need a temporary copy of R that we can decrement.
//
  for ( i = 0; i < m; i++ )
  {
    r2[i] = r[i];
  }

  for ( j = n-1; 0 <= j; j-- )
  {
    i = i4vec_maxloc_last ( m, r2 );

    for ( k = 1; k <= c[j]; k++ )
    {
//
//  By adding 1 rather than setting A(I,J) to 1, we were able to spot
//  an error where the index was "sticking".
//
      a[i+j*m] = a[i+j*m] + 1;

      r2[i] = r2[i] - 1;

      if ( 0 < i )
      {
        i = i - 1;
      }
//
//  There's a special case you have to watch out for.
//  If I was 1, and when you decrement R2(1), I is going to be 1 again,
//  and you're staying in the same column, that's not good.
//
//  The syntax "R2+1" means the vector starting with the second element of R2.
//
      else
      {
        i = i4vec_maxloc_last ( m, r2 );
        if ( i == 0 && k < c[j] )
        {
          i = 1 + i4vec_maxloc_last ( m-1, r2 + 1 );
        }
      }
    }
  }
  delete [] r2;

  return;
}
//****************************************************************************80

void i4mat_01_rowcolsum2 ( int m, int n, int r[], int c[], int a[], 
  bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_01_ROWCOLSUM2 creates a 0/1 I4MAT with given row and column sums.
//
//  Discussion:
//
//    This routine uses network flow optimization to compute the results.
//
//    Given an M vector R and N vector C, there may exist one or more
//    M by N matrices with entries that are 0 or 1, whose row sums are R
//    and column sums are C.
//
//    For convenience, this routine requires that the entries of R and C
//    be given in nonincreasing order.
//
//    There are several requirements on R and C.  The simple requirements
//    are that the entries of R and C must be nonnegative, that the entries
//    of R must each no greater than N, and those of C no greater than M,
//    and that the sum of the entries of R must equal the sum of the 
//    entries of C.
//
//    The final technical requirement is that if we form R*, the conjugate
//    partition of R, then C is majorized by R*, that is, that every partial
//    sum from 1 to K of the entries of C is no bigger than the sum of the same
//    entries of R*, for every K from 1 to N.
//
//    Given these conditions on R and C, there is at least one 0/1 matrix
//    with the given row and column sums.
//
//    The conjugate partition of R is constructed as follows:
//      R*(1) is the number of entries of R that are 1 or greater.
//      R*(2) is the number of entries of R that are 2 or greater.
//      ...
//      R*(N) is the number of entries of R that are N (can't be greater).
//
//  Example:
//
//    M = N = 5
//    R = ( 3, 2, 2, 1, 1 )
//    C = ( 2, 2, 2, 2, 1 )
//
//    A =
//      1 0 1 0 1
//      1 0 0 1 0
//      0 1 0 1 0
//      0 1 0 0 0
//      0 0 1 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//    Jack van Lint, Richard Wilson,
//    A Course in Combinatorics,
//    Oxford, 1992, pages 148-156.
//
//    James Sandeson,
//    Testing Ecobool Patterns,
//    American Scientist,
//    Volume 88, July-August 2000, pages 332-339.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//    These values do not have to be equal.
//
//    Input, int R[M], C[N], the row and column sums desired for the array.
//    Both vectors must be arranged in descending order.
//    The elements of R must be between 0 and N.
//    The elements of C must be between 0 and M.
//    One of the conditions for a solution to exist is that the sum of the
//    elements in R equal the sum of the elements in C.
//
//    Output, int A[M*N], the matrix with the given row and column sums.
//    Each entry of A is 0 or 1.
//
//    Output, bool &ERROR, is true if an error occurred.
//
{
  int *capflo;
  int i;
  int *icut;
  int *iendpt;
  int isink;
  int j;
  int k;
  int nedge;
  int nnode;
  int *node_flow;
  int source;

  error = false;

  capflo = new int[2*2*(m+m*n+n)];
  icut = new int[m+n+2];
  iendpt = new int[2*2*(m+m*n+n)];
  node_flow = new int[m+n+2];
//
//  There are M + N + 2 nodes.  The last two are the special source and sink.
//
  source = m + n + 1;
  isink = m + n + 2;
  nnode = m + n + 2;
//
//  The source is connected to each of the R nodes.
//
  k = 0;

  for ( i = 0; i < m; i++ )
  {
    iendpt[0+2*k] = source;
    iendpt[1+2*k] = i+1;
    capflo[0+2*k] = r[i];
    capflo[1+2*k] = 0;
    k = k + 1;

    iendpt[0+2*k] = i+1;
    iendpt[1+2*k] = source;
    capflo[0+2*k] = r[i];
    capflo[1+2*k] = 0;
    k = k + 1;
  }
//
//  Every R node is connected to every C node, with capacity 1.
//
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      iendpt[0+2*k] = i+1;
      iendpt[1+2*k] = j+1+m;
      capflo[0+2*k] = 1;
      capflo[1+2*k] = 0;
      k = k + 1;

      iendpt[0+2*k] = j+1+m;
      iendpt[1+2*k] = i+1;
      capflo[0+2*k] = 1;
      capflo[1+2*k] = 0;
      k = k + 1;
    }
  }
//
//  Every C node is connected to the sink.
//
  for ( j = 0; j < n; j++ )
  {
    iendpt[0+2*k] = j+1+m;
    iendpt[1+2*k] = isink;
    capflo[0+2*k] = c[j];
    capflo[1+2*k] = 0;
    k = k + 1;

    iendpt[0+2*k] = isink;
    iendpt[1+2*k] = j+1+m;
    capflo[0+2*k] = c[j];
    capflo[1+2*k] = 0;
    k = k + 1;
  }
//
//  Determine the maximum flow on the network.
//
  nedge = k;

  network_flow_max ( nnode, nedge, iendpt, capflo, source, isink, 
    icut, node_flow );
//
//  We have a perfect solution if, and only if, the edges leading from the
//  source, and the edges leading to the sink, are all saturated.
//
  for ( k = 0; k < nedge; k++ )
  {
    i = iendpt[0+2*k];
    j = iendpt[1+2*k] - m;

    if ( i <= m && 1 <= j && j <= n )
    {
      if ( capflo[1+2*k] != 0 && capflo[1+2*k] != 1 )
      {
        error = true;
      }
    }

    if ( iendpt[0+2*k] == source )
    {
      if ( capflo[0+2*k] != capflo[1+2*k] )
      {
        error = true;
      }
    }

    if ( iendpt[1+2*k] == isink )
    {
      if ( capflo[0+2*k] != capflo[1+2*k] )
      {
        error = true;
      }
    }

  }
//
//  If we have a solution, then A(I,J) = the flow on the edge from
//  R node I to C node J.
//
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i+j*m] = 0;
    }
  }

  for ( k = 0; k < nedge; k++ )
  {
    i = iendpt[0+2*k];
    j = iendpt[1+2*k] - m;

    if ( i <= m && 1 <= j && j <= n )
    {
      a[i+j*m] = capflo[1+2*k];
    }
  }

  delete [] icut;
  delete [] iendpt;
  delete [] node_flow;

  return;
}
//****************************************************************************80

void i4mat_perm ( int n, int a[], int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PERM permutes the rows and columns of a square I4MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, int A[N*N].
//    On input, the matrix to be permuted.
//    On output, the permuted matrix.
//
//    Input, int P[N], the permutation.  P(I) is the new number of row
//    and column I.
//
{
  int i;
  int i1;
  int is;
  int it;
  int j;
  int j1;
  int j2;
  int k;
  int lc;
  int nc;
  int temp;

  perm_cycle ( n, p, is, nc, 1 );

  for ( i = 1; i <= n; i++ )
  {
    i1 = - p[i-1];

    if ( 0 < i1 )
    {
      lc = 0;

      for ( ; ; )
      {
        i1 = p[i1-1];
        lc = lc + 1;

        if ( i1 <= 0 )
        {
          break;
        }
      }
      i1 = i;

      for ( j = 1; j <= n; j++ )
      {
        if ( p[j-1] <= 0 )
        {
          j2 = j;
          k = lc;
          for ( ; ; )
          {
            j1 = j2;
            it = a[i1-1+(j1-1)*n];

            for ( ; ; )
            {
              i1 = abs ( p[i1-1] );
              j1 = abs ( p[j1-1] );

              temp = a[i1-1+(j1-1)*n];
              a[i1-1+(j1-1)*n] = it;
              it = temp;

              if ( j1 != j2 )
              {
                continue;
              }

              k = k - 1;

              if ( i1 == i )
              {
                break;
              }
            }

            j2 = abs ( p[j2-1] );

            if ( k == 0 )
            {
              break;
            }
          }
        }
      }
    }
  }
//
//  Restore the positive signs of the data.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = abs ( p[i] );
  }

  return;
}
//****************************************************************************80

void i4mat_perm2 ( int m, int n, int a[], int p[], int q[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PERM2 permutes the rows and columns of a rectangular I4MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int M, number of rows in the matrix.
//
//    Input, int N, number of columns in the matrix.
//
//    Input/output, int A[M*N].
//    On input, the matrix to be permuted.
//    On output, the permuted matrix.
//
//    Input, int P[M], the row permutation.  P(I) is the new number of row I.
//
//    Input, int Q[N].  The column permutation.  Q(I) is the new number
//    of column I.  Note that this routine allows you to pass a single array
//    as both P and Q.
//
{
  int i;
  int i1;
  int is;
  int it;
  int j;
  int j1;
  int j2;
  int k;
  int lc;
  int nc;
  int temp;

  perm_cycle ( m, p, is, nc, 1 );

  if ( 0 < q[0] )
  {
    perm_cycle ( n, q, is, nc, 1 );
  }

  for ( i = 1; i <= m; i++ )
  {
    i1 = -p[i-1];

    if ( 0 < i1 )
    {
      lc = 0;

      for ( ; ; )
      {
        i1 = p[i1-1];
        lc = lc + 1;

        if ( i1 <= 0 )
        {
          break;
        }
      }

      i1 = i;

      for ( j = 1; j <= n; j++ )
      {
        if ( q[j-1] <= 0 )
        {
          j2 = j;
          k = lc;

          for ( ; ; )
          {
            j1 = j2;
            it = a[i1-1+(j1-1)*m];

            for ( ; ; )
            {
              i1 = abs ( p[i1-1] );
              j1 = abs ( q[j1-1] );

              temp = it;
              it = a[i1-1+(j1-1)*m];
              a[i1-1+(j1-1)*m] = temp;

              if ( j1 != j2 )
              {
                continue;
              }

              k = k - 1;

              if ( i1 == i )
              {
                break;
              }
            }

            j2 = abs ( q[j2-1] );

            if ( k == 0 )
            {
              break;
            }
          }
        }
      }
    }
  }
//
//  Restore the positive signs of the data.
//
  for ( i = 0; i < m; i++ )
  {
    p[i] = abs ( p[i] );
  }

  if ( q[0] <= 0 )
  {
    for ( j = 0; j < n; j++ )
    {
      q[j] = abs ( q[j] );
    }
  }
  return;
}
//****************************************************************************80

void i4mat_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT prints an I4MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT_SOME prints some of an I4MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to INCX) entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void i4mat_u1_inverse ( int n, int a[], int b[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_U1_INVERSE inverts a unit upper triangular I4MAT.
//
//  Discussion:
//
//    A unit upper triangular matrix is a matrix with only 1's on the main
//    diagonal, and only 0's below the main diagonal.  Above the main
//    diagonal, the entries may be assigned any value.
//
//    It may be surprising to note that the inverse of an integer unit upper
//    triangular matrix is also an integer unit upper triangular matrix.
//
//    Note that this routine can invert a matrix in place, that is, with no
//    extra storage.  If the matrix is stored in A, then the call
//
//      i4mat_u1_inverse ( n, a, a )
//
//    will result in A being overwritten by its inverse, which can
//    save storage if the original value of A is not needed later.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of rows and columns in the matrix.
//
//    Input, int A[N*N], the unit upper triangular matrix
//    to be inverted.
//
//    Output, int B[N*N], the inverse matrix.
//
{
  int i;
  int isum;
  int j;
  int k;

  for ( j = n; 1 <= j; j-- )
  {
    for ( i = n; 1 <= i; i-- )
    {
      if ( i == j )
      {
        isum = 1;
      }
      else
      {
        isum = 0;
      }

      for ( k = i+1; k <= j; k++ )
      {
        isum = isum - a[i-1+(k-1)*n] * b[k-1+(j-1)*n];
      }
      b[i-1+(j-1)*n] = isum;
    }
  }

  return;
}
//****************************************************************************80

void i4poly ( int n, int a[], int x0, int iopt, int &val )

//****************************************************************************80
//
//  Purpose:
//
//    I4POLY performs operations on I4POLY's in power or factorial form.
//
//  Discussion:
//
//    The power sum form of a polynomial is
//
//      P(X) = A1 + A2*X + A3*X**2 + ... + (AN+1)*X**N
//
//    The Taylor expansion at C has the form
//
//      P(X) = A1 + A2*(X-C) + A3*(X-C)**2 + ... + (AN+1)*(X-C)**N
//
//    The factorial form of a polynomial is
//
//      P(X) = A1 + A2*X + A3*(X)*(X-1) + A4*(X)*(X-1)*(X-2)+...
//        + (AN+1)*(X)*(X-1)*...*(X-N+1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of coefficients in the polynomial
//    (in other words, the polynomial degree + 1)
//
//    Input/output, int A[N], the coefficients of the polynomial.  Depending
//    on the option chosen, these coefficients may be overwritten by those
//    of a different form of the polynomial.
//
//    Input, int X0, for IOPT = -1, 0, or positive, the value of the
//    argument at which the polynomial is to be evaluated, or the
//    Taylor expansion is to be carried out.
//
//    Input, int IOPT, a flag describing which algorithm is to
//    be carried out:
//    -3: Reverse Stirling.  Input the coefficients of the polynomial in
//    factorial form, output them in power sum form.
//    -2: Stirling.  Input the coefficients in power sum form, output them
//    in factorial form.
//    -1: Evaluate a polynomial which has been input in factorial form.
//    0:  Evaluate a polynomial input in power sum form.
//    1 or more:  Given the coefficients of a polynomial in
//    power sum form, compute the first IOPT coefficients of
//    the polynomial in Taylor expansion form.
//
//    Output, int &VAL, for IOPT = -1 or 0, the value of the
//    polynomial at the point X0.
//
{
  int eps;
  int i;
  int m;
  int n1;
  int w;
  int z;

  n1 = i4_min ( n, iopt );
  n1 = i4_max ( 1, n1 );

  if ( iopt < -1 )
  {
    n1 = n;
  }

  eps = i4_max ( -iopt, 0 ) % 2;

  w = - n * eps;

  if ( -2 < iopt )
  {
    w = w + x0;
  }

  for ( m = 1; m <= n1; m++ )
  {
    val = 0;
    z = w;

    for ( i = m; i <= n; i++ )
    {
      z = z + eps;
      val = a[n+m-i-1] + z * val;
      if ( iopt != 0 && iopt != -1 )
      {
        a[n+m-i-1] = val;
      }
    }

    if ( iopt < 0 )
    {
      w = w + 1;
    }
  }

  return;
}
//****************************************************************************80

void i4poly_cyclo ( int n, int phi[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4POLY_CYCLO computes a cyclotomic I4POLY.
//
//  Discussion:
//
//    For 1 <= N, let
//
//      I = SQRT ( - 1 )
//      L = EXP ( 2 * PI * I / N )
//
//    Then the N-th cyclotomic polynomial is defined by
//
//      PHI(N;X) = Product ( 1 <= K <= N and GCD(K,N) = 1 ) ( X - L**K )
//
//    We can use the Moebius MU function to write
//
//      PHI(N;X) = Product ( mod ( D, N ) = 0 ) ( X**D - 1 )**MU(N/D)
//
//    There is a sort of inversion formula:
//
//      X**N - 1 = Product ( mod ( D, N ) = 0 ) PHI(D;X)
//
//  Example:
//
//     N  PHI
//
//     0  1
//     1  X - 1
//     2  X + 1
//     3  X**2 + X + 1
//     4  X**2 + 1
//     5  X**4 + X**3 + X**2 + X + 1
//     6  X**2 - X + 1
//     7  X**6 + X**5 + X**4 + X**3 + X**2 + X + 1
//     8  X**4 + 1
//     9  X**6 + X**3 + 1
//    10  X**4 - X**3 + X**2 - X + 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Raymond Seroul,
//    Programming for Mathematicians,
//    Springer Verlag, 2000, page 269.
//
//  Parameters:
//
//    Input, int N, the index of the cyclotomic polynomial desired.
//
//    Output, int PHI[N+1], the N-th cyclotomic polynomial.
//
{
# define POLY_MAX 100

  int d;
  int den[POLY_MAX+1];
  int den_n;
  int *factor;
  int i;
  int j;
  int mu;
  int nq;
  int nr;
  int num[POLY_MAX+1];
  int num_n;
  int *rem;

  factor = new int[n+1];
  rem = new int[n+1];

  num[0] = 1;
  for ( i = 1; i <= POLY_MAX; i++ )
  {
    num[i] = 0;
  }
  num_n = 0;

  den[0] = 1;
  for ( i = 1; i <= POLY_MAX; i++ )
  {
    den[i] = 0;
  }
  den_n = 0;

  for ( i = 0; i <= n; i++ )
  {
    phi[i] = 0;
  }

  for ( d = 1; d <= n; d++ )
  {
//
//  For each divisor D of N, ...
//
    if ( ( n % d ) == 0 )
    {
      mu = i4_moebius ( n / d );
//
//  ...multiply the numerator or denominator by (X^D-1).
//
      factor[0] = -1;
      for ( j = 1; j <= d-1; j++ )
      {
        factor[j] = 0;
      }
      factor[d] = 1;

      if ( mu == +1 )
      {
        if ( POLY_MAX < num_n + d )
        {
          cerr << "\n";
          cerr << "I4POLY_CYCLO - Fatal error!\n";
          cerr << "  Numerator polynomial degree too high.\n";
          exit ( 1 );
        }

        i4poly_mul ( num_n, num, d, factor, num );

        num_n = num_n + d;
      }
      else if ( mu == -1 )
      {
        if ( POLY_MAX < den_n + d )
        {
          cerr << "\n";
          cerr << "I4POLY_CYCLO - Fatal error!\n";
          cerr << "  Denominator polynomial degree too high.\n";
          exit ( 1 );
        }

        i4poly_mul ( den_n, den, d, factor, den );

        den_n = den_n + d;
      }
    }
  }
//
//  PHI = NUM / DEN
//
  i4poly_div ( num_n, num, den_n, den, nq, phi, nr, rem );

  delete [] factor;
  delete [] rem;

  return;
# undef POLY_MAX
}
//****************************************************************************80

int i4poly_degree ( int na, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4POLY_DEGREE returns the degree of an I4POLY.
//
//  Discussion:
//
//    The degree of a polynomial is the index of the highest power
//    of X with a nonzero coefficient.
//
//    The degree of a constant polynomial is 0.  The degree of the
//    zero polynomial is debatable, but this routine returns the
//    degree as 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the dimension of A.
//
//    Input, int A[NA+1], the coefficients of the polynomials.
//
//    Output, int I4POLY_DEGREE, the degree of the polynomial.
//
{
  int degree;

  degree = na;

  while ( 0 < degree )
  {
    if ( a[degree] != 0 )
    {
      return degree;
    }
    degree = degree - 1;
  }

  return degree;
}
//****************************************************************************80

void i4poly_div ( int na, int a[], int nb, int b[], int &nq, int q[], 
  int &nr, int r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4POLY_DIV computes the quotient and remainder of two I4POLY's.
//
//  Discussion:
//
//    Normally, the quotient and remainder would have rational coefficients.
//    This routine assumes that the special case applies that the quotient
//    and remainder are known beforehand to be integral.
//
//    The polynomials are assumed to be stored in power sum form.
//
//    The power sum form is:
//
//      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the degree of polynomial A.
//
//    Input, int A[NA+1], the coefficients of the polynomial to be divided.
//
//    Input, int NB, the degree of polynomial B.
//
//    Input, int B[NB+1], the coefficients of the divisor polynomial.
//
//    Output, int &NQ, the degree of polynomial Q.
//    If the divisor polynomial is zero, NQ is returned as -1.
//
//    Output, int Q[NA-NB+1], contains the quotient of A/B.
//    If A and B have full degree, Q should be dimensioned Q(0:NA-NB).
//    In any case, Q(0:NA) should be enough.
//
//    Output, int &NR, the degree of polynomial R.
//    If the divisor polynomial is zero, NR is returned as -1.
//
//    Output, int R[NB], contains the remainder of A/B.
//    If B has full degree, R should be dimensioned R(0:NB-1).
//    Otherwise, R will actually require less space.
//
{
  int *a2;
  int i;
  int j;
  int na2;
  int nb2;

  na2 = i4poly_degree ( na, a );

  nb2 = i4poly_degree ( nb, b );

  if ( b[nb2] == 0 )
  {
    nq = -1;
    nr = -1;
    return;
  }

  a2 = new int[na+1];

  for ( i = 0; i <= na2; i++ )
  {
    a2[i] = a[i];
  }

  nq = na2 - nb2;
  nr = nb2 - 1;

  for ( i = nq; 0 <= i; i-- )
  {
    q[i] = a2[i+nb2] / b[nb2];
    a2[i+nb2] = 0;
    for ( j = 0; j < nb2; j++ )
    {
      a2[i+j] = a2[i+j] - q[i] * b[j];
    }
  }

  for ( i = 0; i <= nr; i++ )
  {
    r[i] = a2[i];
  }

  delete [] a2;

  return;
}
//****************************************************************************80

void i4poly_mul ( int na, int a[], int nb, int b[], int c[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4POLY_MUL computes the product of two I4POLY's.
//
//  Discussion:
//
//    The polynomials are in power sum form.
//
//    The power sum form is:
//
//      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the degree of polynomial A.
//
//    Input, int A[NA+1], the coefficients of the first polynomial factor.
//
//    Input, int NB, the degree of polynomial B.
//
//    Input, int B[NB+1], the coefficients of the second polynomial factor.
//
//    Output, int C[NA+NB+1], the coefficients of A * B.
//
{
  int *d;
  int i;
  int j;

  d = new int[na+nb+1];

  for ( i = 0; i <= na+nb; i++ )
  {
    d[i] = 0;
  }

  for ( i = 0; i <= na; i++ )
  {
    for ( j = 0; j <= nb; j++ )
    {
      d[i+j] = d[i+j] + a[i] * b[j];
    }
  }

  for ( i = 0; i <= na+nb; i++ )
  {
    c[i] = d[i];
  }

  delete [] d;

  return;
}
//****************************************************************************80

void i4poly_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4POLY_PRINT prints out an I4POLY.
//
//  Discussion:
//
//    The power sum form is:
//
//      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the degree of polynomial A.
//
//    Input, int A[N+1], the polynomial coefficients.
//    A(0) is the constant term and
//    A(N) is the coefficient of X**N.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int mag;
  int n2;
  char plus_minus;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  n2 = i4poly_degree ( n, a );

  if ( a[n2] < 0 )
  {
    plus_minus = '-';
  }
  else
  {
    plus_minus = ' ';
  }

  mag = abs ( a[n2] );

  if ( 2 <= n2 )
  {
    cout << "p(x) = " << plus_minus << mag << " * x^" << n2 << "\n";
  }
  else if ( n2 == 1 )
  {
    cout << "p(x) = " << plus_minus << mag << " * x" << "\n";
  }
  else if ( n2 == 0 )
  {
    cout << "p(x) = " << plus_minus << mag << "\n";
  }

  for ( i = n2-1; 0 <= i; i-- )
  {
    if ( a[i] < 0.0 )
    {
      plus_minus = '-';
    }
    else
    {
      plus_minus = '+';
    }

    mag = abs ( a[i] );

    if ( mag != 0 )
    {
      if ( 2 <= i )
      {
        cout << "       " << plus_minus << mag << " * x^" << i << "\n";
      }
      else if ( i == 1 )
      {
        cout << "       " << plus_minus << mag << " * x" << "\n";
      }
      else if ( i == 0 )
      {
        cout << "       " << plus_minus << mag << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

int i4poly_to_i4 ( int n, int a[], int x )

//****************************************************************************80
//
//  Purpose:
//
//    I4POLY_TO_I4 evaluates an I4POLY.
//
//  Discussion:
//
//    The power sum form is:
//
//      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the degree of the polynomial.
//
//    Input, int A[N+1], the polynomial coefficients.
//    A[0] is the constant term and
//    A[N] is the coefficient of X**N.
//
//    Input, int X, the point at which the polynomial is to be evaluated.
//
//    Output, int I4POLY_TO_I4, the value of the polynomial.
//
{
  int i;
  int value;

  value = 0;

  for ( i = n; 0 <= i; i-- )
  {
    value = value * x + a[i];
  }

  return value;
}
//****************************************************************************80

bool i4vec_ascends ( int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ASCENDS determines if an I4VEC is (weakly) ascending.
//
//  Example:
//
//    X = ( -8, 1, 2, 3, 7, 7, 9 )
//
//    I4VEC_ASCENDS = TRUE
//
//    The sequence is not required to be strictly ascending.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the array.
//
//    Input, int X[N], the array to be examined.
//
//    Output, bool I4VEC_ASCENDS, is TRUE if the entries of X ascend.
//
{
  int i;

  for ( i = 1; i <= n-1; i++ )
  {
    if ( x[i] < x[i-1] )
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

void i4vec_backtrack ( int n, int maxstack, int stack[], int x[], int &indx, 
  int &k, int &nstack, int ncan[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_BACKTRACK supervises a backtrack search for an I4VEC.
//
//  Discussion:
//
//    The routine tries to construct an integer vector one index at a time,
//    using possible candidates as supplied by the user.
//
//    At any time, the partially constructed vector may be discovered to be
//    unsatisfactory, but the routine records information about where the
//    last arbitrary choice was made, so that the search can be
//    carried out efficiently, rather than starting out all over again.
//
//    First, call the routine with INDX = 0 so it can initialize itself.
//
//    Now, on each return from the routine, if INDX is:
//      1, you've just been handed a complete candidate vector;
//         Admire it, analyze it, do what you like.
//      2, please determine suitable candidates for position X(K).
//         Return the number of candidates in NCAN(K), adding each
//         candidate to the end of STACK, and increasing NSTACK.
//      3, you're done.  Stop calling the routine;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 July 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of positions to be filled in the vector.
//
//    Input, int MAXSTACK, the maximum length of the stack.
//
//    Input, int STACK[MAXSTACK], a list of all current candidates for
//    all positions 1 through K.
//
//    Input/output, int X[N], the partial or complete candidate vector.
//
//    Input/output, int &INDX, a communication flag.
//    On input,
//      0 to start a search.
//    On output:
//      1, a complete output vector has been determined and returned in X(1:N);
//      2, candidates are needed for position X(K);
//      3, no more possible vectors exist.
//
//    Input/output, int &K, if INDX=2, the current vector index being considered.
//
//    Input/output, int &NSTACK, the current length of the stack.
//
//    Input/output, int NCAN[N], lists the current number of candidates for
//    positions 1 through K.
//
{
//
//  If this is the first call, request a candidate for position 1.
//
  if ( indx == 0 )
  {
    k = 1;
    nstack = 0;
    indx = 2;
    return;
  }
//
//  Examine the stack.
//
  for ( ; ; )
  {
//
//  If there are candidates for position K, take the first available
//  one off the stack, and increment K.
//
//  This may cause K to reach the desired value of N, in which case
//  we need to signal the user that a complete set of candidates
//  is being returned.
//
    if ( 0 < ncan[k-1] )
    {
      x[k-1] = stack[(nstack)-1];
      nstack = nstack - 1;

      ncan[k-1] = ncan[k-1] - 1;

      if ( k != n )
      {
        k = k + 1;
        indx = 2;
      }
      else
      {
        indx = 1;
      }
      break;
    }
//
//  If there are no candidates for position K, then decrement K.
//  If K is still positive, repeat the examination of the stack.
//
    else
    {
      k = k - 1;

      if ( k <= 0 )
      {
        indx = 3;
        break;
      }
    }
  }
  return;
}
//****************************************************************************80

bool i4vec_descends ( int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_DESCENDS determines if an I4VEC is decreasing.
//
//  Example:
//
//    X = ( 9, 7, 7, 3, 2, 1, -8 )
//
//    I4VEC_DESCENDS = TRUE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the array.
//
//    Input, int X[N], the array to be examined.
//
//    Output, bool I4VEC_DESCEND, is TRUE if the entries of the array descend.
//
{
  int i;

  for ( i = 0; i < n-1; i++ )
  {
    if ( x[i] < x[i+1] )
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

int i4vec_frac ( int n, int a[], int k )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_FRAC searches for the K-th smallest entry in an I4VEC.
//
//  Discussion:
//
//    Hoare's algorithm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input/output, int A[N].
//    On input, A is the array to search.
//    On output, the elements of A have been somewhat rearranged.
//
//    Input, int K, the fractile to be sought.  If K = 1, the minimum
//    entry is sought.  If K = N, the maximum is sought.  Other values
//    of K search for the entry which is K-th in size.  K must be at
//    least 1, and no greater than N.
//
//    Output, int I4VEC_FRAC, the value of the K-th fractile of A.
//
{
  int afrac;
  int i;
  int iryt;
  int j;
  int left;
  int temp;
  int x;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal nonpositive value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( k <= 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal nonpositive value of K = " << k << "\n";
    exit ( 1 );
  }

  if ( n < k )
  {
    cerr << "\n";
    cerr << "I4VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal N < K, K = " << k << "\n";
    exit ( 1 );
  }

  left = 1;
  iryt = n;

  for ( ; ; )
  {
    if ( iryt <= left )
    {
      afrac = a[k-1];
      break;
    }

    x = a[k-1];
    i = left;
    j = iryt;

    for ( ; ; )
    {
      if ( j < i )
      {
        if ( j < k )
        {
          left = i;
        }
        if ( k < i )
        {
          iryt = j;
        }
        break;
      }
//
//  Find I so that X <= A(I).
//
      while ( a[i-1] < x )
      {
        i = i + 1;
      }
//
//  Find J so that A(J) <= X
//
      while ( x < a[j-1] )
      {
        j = j - 1;
      }

      if ( i <= j )
      {
        temp = a[i-1];
        a[i-1] = a[j-1];
        a[j-1] = temp;
        i = i + 1;
        j = j - 1;
      }
    }
  }
  return afrac;
}
//****************************************************************************80

void i4vec_heap_d ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_HEAP_D reorders an I4VEC into a descending heap.
//
//  Discussion:
//
//    A heap is an array A with the property that, for every index J,
//    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
//    2*J+1 and 2*J+2 are legal).
//
//  Diagram:
//
//                  A(0)
//                /      \
//            A(1)         A(2)
//          /     \        /  \
//      A(3)       A(4)  A(5) A(6)
//      /  \       /   \
//    A(7) A(8)  A(9) A(10)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the input array.
//
//    Input/output, int A[N].
//    On input, an unsorted array.
//    On output, the array has been reordered into a heap.
//
{
  int i;
  int ifree;
  int key;
  int m;
//
//  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
//
  for ( i = (n/2)-1; 0 <= i; i-- )
  { 
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = a[i];
    ifree = i;

    for ( ; ; )
    {
//
//  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
//  IFREE.  (One or both may not exist because they equal or exceed N.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n <= m )
      {
        break;
      }
      else
      {
//
//  Does the second position exist?
//
        if ( m + 1 < n )
        {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }
      }
    }
//
//  When you have stopped shifting items up, return the item you
//  pulled out back to the heap.
//
    a[ifree] = key;
  }

  return;
}
//****************************************************************************80

int i4vec_index ( int n, int a[], int aval )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDEX returns the location of the first occurrence of a given value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be searched.
//
//    Input, int AVAL, the value to be indexed.
//
//    Output, int I4VEC_INDEX, the first location in A which has the
//    value AVAL, or -1 if no such index exists.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] == aval )
    {
      return i;
    }
  }
  return -1;
}
//****************************************************************************80

void i4vec_indicator ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR sets an I4VEC to the indicator vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int A[N], the initialized array.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }
  return;
}
//****************************************************************************80

int *i4vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR_NEW sets an I4VEC to the indicator vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int I4VEC_INDICATOR_NEW[N], the initialized array.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }
  return a;
}
//****************************************************************************80

int i4vec_max ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MAX returns the value of the maximum element in an I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MIN, the value of the maxnimum element.  This
//    is set to 0 if N <= 0.
//
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < a[i] )
    {
      value = a[i];
    }
  }
  return value;
}
//****************************************************************************80

int i4vec_maxloc_last ( int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MAXLOC_LAST returns the index of the last maximal I4VEC entry.
//
//  Example:
//
//    X = ( 5, 1, 2, 5, 0, 5, 3 )
//
//    I4VEC_MAXLOC_LAST = 5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the array.
//
//    Input, int X[N], the array to be examined.
//
//    Output, int I4VEC_MAXLOC_LAST, the index of the last element of
//    X of maximal value.
//
{
  int i;
  int index;
  int value;

  index = 0;
  value = x[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value <= x[i] )
    {
      index = i;
      value = x[i];
    }
  }
  return index;
}
//****************************************************************************80

int i4vec_min ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MIN returns the value of the minimum element in an I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MIN, the value of the minimum element.  This
//    is set to 0 if N <= 0.
//
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] < value )
    {
      value = a[i];
    }
  }
  return value;
}
//****************************************************************************80

bool i4vec_pairwise_prime ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PAIRWISE_PRIME checks whether an I4VEC is pairwise prime.
//
//  Discussion:
//
//    Two positive integers I and J are pairwise prime if they have no common
//    factor greater than 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values to check.
//
//    Input, int A[N], the vector of integers.
//
//    Output, bool I4VEC_PAIRWISE_PRIME, is TRUE if the vector of integers
//    is pairwise prime.
//
{
  int i;
  int j;

  for ( i = 0; i < n; i++ )
  {
    for ( j = i+1; j < n; j++ )
    {
      if ( i4_gcd ( a[i], a[j] ) != 1 ) 
      {
        return false;
      }
    }
  }
  return true;
}
//****************************************************************************80

int i4vec_product ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRODUCT multiplies the entries of an I4VEC.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_PRODUCT = 24
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector
//
//    Output, int I4VEC_PRODUCT, the product of the entries of A.
//
{
  int i;
  int product;

  product = 1;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}
//****************************************************************************80

void i4vec_reverse ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_REVERSE reverses the elements of an integer vector.
//
//  Example:
//
//    Input:
//
//      N = 5,
//      A = ( 11, 12, 13, 14, 15 ).
//
//    Output:
//
//      A = ( 15, 14, 13, 12, 11 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A[N], the array to be reversed.
//
{
  int i;
  int temp;

  for ( i = 0; i < n/2; i++ )
  {
    temp = a[i];
    a[i] = a[n-1-i];
    a[n-1-i] = temp;
  }
  return;
}
//****************************************************************************80

void i4vec_sort_bubble_a ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_BUBBLE_A ascending sorts an integer array using bubble sort.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  int i;
  int j;
  int k;

  for ( i = 0; i < n-1; i++ )
  {
    for ( j = i+1; j < n; j++ )
    {
      if ( a[j] < a[i] )
      {
        k = a[i];
        a[i] = a[j];
        a[j] = k;
      }
    }
  }
  return;
}
//****************************************************************************80

void i4vec_sort_heap_a ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  int n1;
  int temp;

  if ( n <= 1 )
  {
    return;
  }
//
//  1: Put A into descending heap form.
//
  i4vec_heap_d ( n, a );
//
//  2: Sort A.
//
//  The largest object in the heap is in A[0].
//  Move it to position A[N-1].
//
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
//
//  Consider the diminished heap of size N1.
//
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
//
//  Restore the heap structure of the initial N1 entries of A.
//
    i4vec_heap_d ( n1, a );
//
//  Take the largest object from A[0] and move it to A[N1-1].
//
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }
  return;
}
//****************************************************************************80

int *i4vec_sort_heap_index_a ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an integer vector.
//
//  Discussion:
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      A(INDX(I)), I = 1 to N is sorted,
//
//    after which A(I), I = 1 to N is sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], an array to be index-sorted.
//
//    Output, int I4VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  int aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  indx = i4vec_indicator_new ( n );

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt-1];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt-1];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        for ( i = 0; i < n; i++ )
        {
          indx[i] = indx[i] - 1;
        }
        return indx;
      }
    }
    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]-1] < a[indx[j]-1] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]-1] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
}
//****************************************************************************80

int *i4vec_sort_heap_index_d ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an I4VEC.
//
//  Discussion:
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      A(INDX(I)), I = 1 to N is sorted,
//
//    after which A(I), I = 1 to N is sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], an array to be index-sorted.
//
//    Output, int I4VEC_SORT_HEAP_INDEX_D[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  int aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  indx = i4vec_indicator_new ( n );

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt-1];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt-1];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        for ( i = 0; i < n; i++ )
        {
          indx[i] = indx[i] - 1;
        }
        return indx;
      }

    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j]-1] < a[indx[j-1]-1] )
        {
          j = j + 1;
        }
      }

      if ( a[indx[j-1]-1] < aval )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }

    }

    indx[i-1] = indxt;

  }
}
//****************************************************************************80

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

int *i4vec_uniform ( int n, int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
//
//  Discussion:
//
//    The pseudorandom numbers should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, integer N, the dimension of the vector.
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int IVEC_UNIFORM[N], a vector of random values between A and B.
//
{
  int i;
  int k;
  float r;
  int value;
  int *x;
  
  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  x = new int[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + 2147483647;
    }

    r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
    r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 ) 
      +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
    value = r4_nint ( r );

    value = i4_max ( value, i4_min ( a, b ) );
    value = i4_min ( value, i4_max ( a, b ) );

    x[i] = value;
  }

  return x;
}
//****************************************************************************80

void i4vec0_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC0_PRINT prints an integer vector (0-based indices).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n-1; i++ ) 
  {
    cout << setw(6) << i    << "  " 
         << setw(8) << a[i] << "\n";
  }

  return;
}
//****************************************************************************80

void i4vec1_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC1_PRINT prints an integer vector (one-based indices).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n-1; i++ ) 
  {
    cout                    << "  "
         << setw(6) << i+1  << "  " 
         << setw(8) << a[i] << "\n";
  }

  return;
}
//****************************************************************************80

void index_box2_next_2d ( int n1, int n2, int ic, int jc, int &i, int &j, 
  bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_BOX2_NEXT_2D produces index vectors on the surface of a box in 2D.
//
//  Discussion:
//
//    The box has center at (IC,JC), and half-widths N1 and N2.
//    The index vectors are exactly those which are between (IC-N1,JC-N1) and
//    (IC+N1,JC+N2) with the property that at least one of I and J
//    is an "extreme" value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the half-widths of the box, that is, the
//    maximum distance allowed between (IC,JC) and (I,J).
//
//    Input, int IC, JC, the central cell of the box.
//
//    Input/output, int &I, &J.  On input, the previous index set.
//    On output, the next index set.  On the first call, MORE should
//    be set to FALSE, and the input values of I and J are ignored.
//
//    Input/output, bool &MORE.
//    On the first call for a given box, the user should set MORE to FALSE.
//    On return, the routine sets MORE to TRUE.
//    When there are no more indices, the routine sets MORE to FALSE.
//
{
  if ( !( more ) )
  {
    more = true;
    i = ic - n1;
    j = jc - n2;
    return;
  }

  if ( i == ic + n1 && j == jc + n2 )
  {
    more = false;
    return;
  }
//
//  Increment J.
//
  j = j + 1;
//
//  Check J.
//
  if ( jc + n2 < j )
  {
    j = jc - n2;
    i = i + 1;
  }
  else if ( j < jc + n2 && ( i == ic - n1 || i == ic + n1 ) )
  {
    return;
  }
  else
  {
    j = jc + n2;
    return;
  }

  return;
}
//****************************************************************************80

void index_box2_next_3d ( int n1, int n2, int n3, int ic, int jc, int kc, 
  int &i, int &j, int &k, bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_BOX2_NEXT_3D produces index vectors on the surface of a box in 3D.
//
//  Discussion:
//
//    The box has a central cell of (IC,JC,KC), with half widths of
//    (N1,N2,N3).  The index vectors are exactly those between
//    (IC-N1,JC-N2,KC-N3) and (IC+N1,JC+N2,KC+N3) with the property that 
//    at least one of I, J, and K is an "extreme" value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the "half widths" of the box, that is, the
//    maximum distances from the central cell allowed for I, J and K.
//
//    Input, int IC, JC, KC, the central cell of the box.
//
//    Input/output, int &I, &J, &K.  On input, the previous index set.
//    On output, the next index set.  On the first call, MORE should
//    be set to FALSE, and the input values of I, J, and K are ignored.
//
//    Input/output, bool &MORE.
//    On the first call for a given box, the user should set MORE to FALSE.
//    On return, the routine sets MORE to TRUE.
//    When there are no more indices, the routine sets MORE to FALSE.
//
{
  if ( !more )
  {
    more = true;
    i = ic - n1;
    j = jc - n2;
    k = kc - n3;
    return;
  }

  if ( i == ic + n1 && j == jc + n2 && k == kc + n3 ) 
  {
    more = false;
    return;
  }
//
//  Increment K.
//
  k = k + 1;
//
//  Check K.
//
  if ( kc + n3 < k )
  {
    k = kc - n3;
    j = j + 1;
  }
  else if ( k < kc + n3 &&
    ( i == ic - n1 || i == ic + n1 ||
      j == jc - n2 || j == jc + n2 ) )
  {
    return;
  }
  else
  {
    k = kc + n3;
    return;
  }
//
//  Check J.
//
  if ( jc + n2 < j )
  {
    j = jc - n2;
    i = i + 1;
  }
  else if ( j < jc + n2 &&
    ( i == ic - n1 || i == ic + n1 || 
      k == kc - n3 || k == kc + n3 ) )
  {
    return;
  }
  else
  {
    j = jc + n2;
    return;
  }

  return;
}
//****************************************************************************80

void index_box_next_2d ( int n1, int n2, int &i, int &j, bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_BOX_NEXT_2D produces index vectors on the surface of a box in 2D.
//
//  Discussion:
//
//    The index vectors are exactly those which are between (1,1) and
//    (N1,N2) with the property that at least one of I and J
//    is an "extreme" value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the "dimensions" of the box, that is, the
//    maximum values allowed for I and J.  The minimum values are
//    assumed to be 1.
//
//    Input/output, int &I, &J.  On input, the previous index set.
//    On output, the next index set.  On the first call, MORE should
//    be set to FALSE, and the input values of I and J are ignored.
//
//    Input/output, bool &MORE.
//    On the first call for a given box, the user should set MORE to FALSE.
//    On return, the routine sets MORE to TRUE.
//    When there are no more indices, the routine sets MORE to FALSE.
//
{
  if ( !more )
  {
    more = true;
    i = 1;
    j = 1;
    return;
  }

  if ( i == n1 && j == n2 )
  {
    more = false;
    return;
  }
//
//  Increment J.
//
  j = j + 1;
//
//  Check J.
//
  if ( n2 < j )
  {
    j = 1;
    i = i + 1;
  }
  else if ( j < n2 && ( i == 1 || i == n1 ) )
  {
    return;
  }
  else
  {
    j = n2;
    return;
  }

  return;
}
//****************************************************************************80

void index_box_next_3d ( int n1, int n2, int n3, int &i, int &j, int &k, 
  bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_BOX_NEXT_3D produces index vectors on the surface of a box in 3D.
//
//  Discussion:
//
//    The index vectors are exactly those which are between (1,1,1) and
//    (N1,N2,N3) with the property that at least one of I, J, and K
//    is an "extreme" value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the "dimensions" of the box, that is, the
//    maximum values allowed for I, J and K.  The minimum values are
//    assumed to be 1.
//
//    Input/output, int &I, &J, &K.  On input, the previous index set.
//    On output, the next index set.  On the first call, MORE should
//    be set to FALSE, and the input values of I, J, and K are ignored.
//
//    Input/output, bool &MORE.
//    On the first call for a given box, the user should set MORE to FALSE.
//    On return, the routine sets MORE to TRUE.
//    When there are no more indices, the routine sets MORE to FALSE.
//
{
  if ( !more )
  {
    more = true;
    i = 1;
    j = 1;
    k = 1;
    return;
  }

  if ( i == n1 && j == n2 && k == n3 )
  {
    more = false;
    return;
  }
//
//  Increment K.
//
  k = k + 1;
//
//  Check K.
//
  if ( n3 < k )
  {
    k = 1;
    j = j + 1;
  }
  else if ( k < n3 && ( i == 1 || i == n1 || j == 1 || j == n2 ) )
  {
    return;
  }
  else
  {
    k = n3;
    return;
  }
//
//  Check J.
//
  if ( n2 < j )
  {
    j = 1;
    i = i + 1;
  }
  else if ( j < n2 && ( i == 1 || i == n1 || k == 1 || k == n3 ) )
  {
    return;
  }
  else
  {
    j = n2;
    return;
  }

  return;
}
//****************************************************************************80

void index_next0 ( int n, int hi, int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_NEXT0 generates all index vectors within given upper limits.
//
//  Discussion:
//
//    The index vectors are generated in such a way that the reversed
//    sequences are produced in lexicographic order.
//
//  Example:
//
//    N = 3,
//    HI = 3
//
//    1   2   3
//    ---------
//    1   1   1
//    2   1   1
//    3   1   1
//    1   2   1
//    2   2   1
//    3   2   1
//    1   3   1
//    2   3   1
//    3   3   1
//    1   1   2
//    2   1   2
//    3   1   2
//    1   2   2
//    2   2   2
//    3   2   2
//    1   3   2
//    2   3   2
//    3   3   2
//    1   1   3
//    2   1   3
//    3   1   3
//    1   2   3
//    2   2   3
//    3   2   3
//    1   3   3
//    2   3   3
//    3   3   3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, int HI, the upper limit for the array indices.
//    The lower limit is implicitly 1 and HI must be at least 1.
//
//    Input/output, int A[N].
//    On startup calls, with MORE = FALSE, the input value of A
//    doesn't matter, because the routine initializes it.
//    On calls with MORE = TRUE, the input value of A must be
//    the output value of A from the previous call.  (In other words,
//    just leave it alone!).
//    On output, A contains the successor set of indices to the input
//    value.
//
//    Input/output, bool &MORE.  Set this variable FALSE before
//    the first call.  Normally, MORE will be returned TRUE but
//    once all the vectors have been generated, MORE will be
//    reset FALSE and you should stop calling the program.
//
{
  int i;
  int inc;

  if ( !more )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 1;
    }

    if ( hi < 1 )
    {
      more = false;
      cerr << "\n";
      cerr << "INDEX_NEXT0 - Fatal error!\n";
      cerr << "  HI is " << hi << "\n";
      cerr << "  but HI must be at least 1.\n";
      exit ( 1 );
    }
  }
  else
  {
    inc = 0;

    while ( hi <= a[inc] )
    {
      a[inc] = 1;
      inc = inc + 1;
    }

    a[inc] = a[inc] + 1;
  }
//
//  See if there are more entries to compute.
//
  more = false;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < hi )
    {
      more = true;
    }
  }

  return;
}
//****************************************************************************80

void index_next1 ( int n, int hi[], int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_NEXT1 generates all index vectors within given upper limits.
//
//  Discussion:
//
//    The index vectors are generated in such a way that the reversed
//    sequences are produced in lexicographic order.
//
//  Example:
//
//    N = 3,
//    HI(1) = 4, HI(2) = 2, HI(3) = 3
//
//    1   2   3
//    ---------
//    1   1   1
//    2   1   1
//    3   1   1
//    4   1   1
//    1   2   1
//    2   2   1
//    3   2   1
//    4   2   1
//    1   1   2
//    2   1   2
//    3   1   2
//    4   1   2
//    1   2   2
//    2   2   2
//    3   2   2
//    4   2   2
//    1   1   3
//    2   1   3
//    3   1   3
//    4   1   3
//    1   2   3
//    2   2   3
//    3   2   3
//    4   2   3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, int HI[N], the upper limits for the array indices.
//    The lower limit is implicitly 1, and each HI(I) should be at least 1.
//
//    Input/output, int A[N].
//    On startup calls, with MORE = FALSE, the input value of A
//    doesn't matter, because the routine initializes it.
//    On calls with MORE = TRUE, the input value of A must be
//    the output value of A from the previous call.  (In other words,
//    just leave it alone!).
//    On output, A contains the successor set of indices to the input
//    value.
//
//    Input/output, bool &MORE.  Set this variable FALSE before
//    the first call.  Normally, MORE will be returned TRUE but
//    once all the vectors have been generated, MORE will be
//    reset FALSE and you should stop calling the program.
//
{
  int i;
  int inc;

  if ( !more )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 1;
    }

    for ( i = 0; i < n; i++ )
    {
      if ( hi[i] < 1 )
      {
        more = false;
        cerr << "\n";
        cerr << "INDEX_NEXT1 - Fatal error!\n";
        cerr << "  Entry " << i << " of HI is " << hi[i] << "\n";
        cerr << "  but all entries must be at least 1.\n";
        exit ( 1 );
      }
    }
  }
  else
  {
    inc = 0;

    while ( hi[inc] <= a[inc] )
    {
      a[inc] = 1;
      inc = inc + 1;
    }

    a[inc] = a[inc] + 1;
  }
//
//  See if there are more entries to compute.
//
  more = false;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < hi[i] )
    {
      more = true;
    }
  }

  return;
}
//****************************************************************************80

void index_next2 ( int n, int lo[], int hi[], int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_NEXT2 generates all index vectors within given lower and upper limits.
//
//  Example:
//
//    N = 3,
//    LO(1) = 1, LO(2) = 10, LO(3) = 4
//    HI(1) = 2, HI(2) = 11, HI(3) = 6
//
//    1   2   3
//    ---------
//    1  10   4
//    2  10   4
//    1  11   4
//    2  11   4
//    1  10   5
//    2  10   5
//    1  11   5
//    2  11   5
//    1  10   6
//    2  10   6
//    1  11   6
//    2  11   6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.  The rank of
//    the object being indexed.
//
//    Input, int LO[N], HI[N], the lower and upper limits for the array
//    indices.  LO(I) should be less than or equal to HI(I), for each I.
//
//    Input/output, int A[N].
//    On startup calls, with MORE = FALSE, the input value of A
//    doesn't matter, because the routine initializes it.
//    On calls with MORE = TRUE, the input value of A must be
//    the output value of A from the previous call.  (In other words,
//    just leave it alone!).
//    On output, A contains the successor set of indices to the input
//    value.
//
//    Input/output, bool &MORE.  Set this variable FALSE before
//    the first call.  Normally, MORE will be returned TRUE but
//    once all the vectors have been generated, MORE will be
//    reset FALSE and you should stop calling the program.
//
{
  int i;
  int inc;

  if ( !more )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = lo[i];
    }

    for ( i = 0; i < n; i++ )
    {
      if ( hi[i] < lo[i] )
      {
        more = false;
        cerr << "\n";
        cerr << "INDEX_NEXT2 - Fatal error!\n";
        cerr << "  Entry " << i << " of HI is " << hi[i] << "\n";
        cerr << "  Entry " << i << " of LO is " << lo[i] << "\n";
        cerr << "  but LO(I) <= HI(I) is required.\n";
        exit ( 1 );
      }
    }
  }
  else
  {
    inc = 0;

    while ( hi[inc] <= a[inc] )
    {
      a[inc] = lo[inc];
      inc = inc + 1;
    }

    a[inc] = a[inc] + 1;
  }
//
//  See if there are more entries to compute.
//
  more = false;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < hi[i] )
    {
      more = true;
    }
  }

  return;
}
//****************************************************************************80

int index_rank0 ( int n, int hi, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_RANK0 ranks an index vector within given upper limits.
//
//  Example:
//
//    N = 3,
//    HI = 3
//    A = ( 3, 1, 2 )
//
//    RANK = 12
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, int HI, the upper limit for the array indices.
//    The lower limit is implicitly 1, and HI should be at least 1.
//
//    Input, int A[N], the index vector to be ranked.
//
//    Output, int INDEX_RANK0, the rank of the index vector, or -1 if A
//    is not a legal index.
//
{
  int i;
  int range;
  int rank;

  rank = -1;
  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < 1 || hi < a[i] )
    {
      return rank;
    }
  }

  rank = 0;
  for ( i = n-1; 0 <= i; i-- )
  {
    rank = hi * rank + a[i];
  }

  rank = 1;
  range = 1;
  for ( i = 0; i < n; i++ )
  {
    rank = rank + ( a[i] - 1 ) * range;
    range = range * hi;
  }

  return rank;
}
//****************************************************************************80

int index_rank1 ( int n, int hi[], int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_RANK1 ranks an index vector within given upper limits.
//
//  Example:
//
//    N = 3,
//    HI(1) = 4, HI(2) = 2, HI(3) = 3
//    A = ( 4, 1, 2 )
//
//    RANK = 12
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, int HI[N], the upper limits for the array indices.
//    The lower limit is implicitly 1, and each HI(I) should be at least 1.
//
//    Input, int A[N], the index to be ranked.
//
//    Output, int INDEX_RANK1, the rank of the index vector, or -1 if A
//    is not a legal index.
//
{
  int i;
  int range;
  int rank;

  rank = -1;
  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < 1 || hi[i] < a[i] )
    {
      return rank;
    }
  }

  rank = 0;
  for ( i = n-1; 0 <= i; i-- )
  {
    rank = hi[i] * rank + a[i];
  }

  rank = 1;
  range = 1;
  for ( i = 0; i < n; i++ )
  {
    rank = rank + ( a[i] - 1 ) * range;
    range = range * hi[i];
  }

  return rank;
}
//****************************************************************************80

int index_rank2 ( int n, int lo[], int hi[], int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_RANK2 ranks an index vector within given lower and upper limits.
//
//  Example:
//
//    N = 3,
//    LO(1) = 1, LO(2) = 10, LO(3) = 4
//    HI(1) = 2, HI(2) = 11, HI(3) = 6
//    A = ( 1, 11, 5 )
//
//    RANK = 7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, int LO[N], HI[N], the lower and upper limits for the array
//    indices.  LO(I) should be less than or equal to HI(I), for each I.
//
//    Input, int A[N], the index vector to be ranked.
//
//    Output, int INDEX_RANK2, the rank of the index vector, or -1 if A
//    is not a legal index vector.
//
{
  int i;
  int range;
  int rank;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < lo[i] || hi[i] < a[i] )
    {
      rank = -1;
      return rank;
    }
  }

  rank = 1;
  range = 1;
  for ( i = 0; i < n; i++ )
  {
    rank = rank + ( a[i] - lo[i] ) * range;
    range = range * ( hi[i] + 1 - lo[i] );
  }

  return rank;
}
//****************************************************************************80

void index_unrank0 ( int n, int hi, int rank, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_UNRANK0 unranks an index vector within given upper limits.
//
//  Example:
//
//    N = 3,
//    HI = 3
//    RANK = 12
//
//    A = ( 3, 1, 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, int HI, the upper limit for the array indices.
//    The lower limit is implicitly 1, and HI should be at least 1.
//
//    Input, int RANK, the rank of the desired index vector.
//
//    Output, int A[N], the index vector of the given rank.
//
{
  int i;
  int j;
  int k;
  int range;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
//
//  The rank might be too small.
//
  if ( rank < 1 )
  {
    return;
  }

  range = ( int ) pow ( ( double ) hi, n );
//
//  The rank might be too large.
//
  if ( range < rank )
  {
    return;
  }

  k = rank - 1;

  for ( i = n-1; 0 <= i; i-- )
  {
    range = range / hi;
    j = k / range;
    a[i] = j + 1;
    k = k - j * range;
  }

  return;
}
//****************************************************************************80

void index_unrank1 ( int n, int hi[], int rank, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_UNRANK1 unranks an index vector within given upper limits.
//
//  Example:
//
//    N = 3,
//    HI(1) = 4, HI(2) = 2, HI(3) = 3
//    RANK = 11
//
//    A = ( 3, 1, 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, int HI[N], the upper limits for the array indices.
//    The lower limit is implicitly 1, and each HI(I) should be at least 1.
//
//    Input, int RANK, the rank of the desired index vector.
//
//    Output, int A[N], the index vector of the given rank.
//
{
  int i;
  int j;
  int k;
  int range;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
//
//  The rank might be too small.
//
  if ( rank < 1 )
  {
    return;
  }

  range = i4vec_product ( n, hi );
//
//  The rank might be too large.
//
  if ( range < rank )
  {
    return;
  }

  k = rank - 1;

  for ( i = n-1; 0 <= i; i-- )
  {
    range = range / hi[i];
    j = k / range;
    a[i] = j + 1;
    k = k - j * range;
  }

  return;
}
//****************************************************************************80

void index_unrank2 ( int n, int lo[], int hi[], int rank, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_UNRANK2 unranks an index vector within given lower and upper limits.
//
//  Example:
//
//    N = 3,
//    LO(1) = 1, LO(2) = 10, LO(3) = 4
//    HI(1) = 2, HI(2) = 11, HI(3) = 6
//    RANK = 7
//
//    A = ( 1, 11, 5 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, int LO[N], HI[N], the lower and upper limits for the array
//    indices.  It should be the case that LO(I) <= HI(I) for each I.
//
//    Input, int RANK, the rank of the desired index.
//
//    Output, int A[N], the index vector of the given rank.
//
{
  int i;
  int j;
  int k;
  int range;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
//
//  The rank might be too small.
//
  if ( rank < 1 )
  {
    return;
  }

  range = 1;
  for ( i = 0; i < n; i++ )
  {
    range = range * ( hi[i] + 1 - lo[i] );
  }
//
//  The rank might be too large.
//
  if ( range < rank )
  {
    return;
  }

  k = rank - 1;
  for ( i = n-1; 0 <= i; i-- )
  {
    range = range / ( hi[i] + 1 - lo[i] );
    j = k / range;
    a[i] = j + lo[i];
    k = k - j * range;
  }

  return;
}
//****************************************************************************80

void ins_perm ( int n, int ins[], int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    INS_PERM computes a permutation from its inversion sequence.
//
//  Discussion:
//
//    For a given permutation P acting on objects 1 through N, the
//    inversion sequence INS is defined as:
//
//      INS(1) = 0
//      INS(I) = number of values J < I for which P(I) < P(J).
//
//  Example:
//
//    Input:
//
//      ( 0, 0, 2, 1, 3 )
//
//    Output:
//
//      ( 3, 5, 1, 4, 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input, int INS[N], the inversion sequence of a permutation.
//    It must be the case that 0 <= INS(I) < I for 1 <= I <= N.
//
//    Output, int P[N], the permutation.
//
{
  int i;
  int itemp;
  int j;

  i4vec_indicator ( n, p );

  for ( i = n-1; 1 <= i; i-- )
  {
    itemp = p[i-ins[i]];

    for ( j = i-ins[i]; j <= i-1; j++ )
    {
      p[j] = p[j+1];
    }

    p[i] = itemp;

  }

  return;
}
//****************************************************************************80

int inverse_mod_n ( int b, int n )

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_MOD_N computes the inverse of B mod N.
//
//  Discussion:
//
//    If 
//
//      Y = inverse_mod_n ( B, N )
//
//    then
//
//      mod ( B * Y, N ) = 1
//
//    The value Y will exist if and only if B and N are relatively prime.
//
//  Examples:
//
//    B  N  Y
//
//    1  2  1
//
//    1  3  1
//    2  3  2
//
//    1  4  1
//    2  4  0
//    3  4  3
//
//    1  5  1
//    2  5  3
//    3  5  2
//    4  5  4
//
//    1  6  1
//    2  6  0
//    3  6  0
//    4  6  0
//    5  6  5
//
//    1  7  1
//    2  7  4
//    3  7  5
//    4  7  2
//    5  7  3
//    6  7  6
//
//    1  8  1
//    2  8  0
//    3  8  3
//    4  8  0
//    5  8  5
//    6  8  0
//    7  8  7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int B, the number whose inverse mod N is desired.
//    B should be positive.  Normally, B < N, but this is not required.
//
//    Input, int N, the number with respect to which the
//    modulus is computed.  N should be positive.
//
//    Output, int INVERSE_MOD_N, the inverse of B mod N, or 0 if there
//    is not inverse for B mode N.
//
{
  int b0;
  int n0;
  int q;
  int r;
  int t;
  int t0;
  int temp;
  int y;

  n0 = n;
  b0 = b;
  t0 = 0;
  t = 1;

  q = n / b;
  r = n - q * b;

  while ( 0 < r )
  {
    temp = t0 - q * t;

    if ( 0 <= temp )
    {
      temp = ( temp % n );
    }

    if ( temp < 0 )
    {
      temp = n - ( ( - temp ) % n );
    }

    t0 = t;
    t = temp;
    n0 = b0;
    b0 = r;
    q = n0 / b0;
    r = n0 - q * b0;
  }

  if ( b0 != 1 )
  {
    y = 0;
    return y;
  }

  y = ( t % n );

  return y;
}
//****************************************************************************80

void involute_enum ( int n, int s[] )

//****************************************************************************80
//
//  Purpose:
//
//    INVOLUTE_ENUM enumerates the involutions of N objects.
//
//  Discussion:
//
//    An involution is a permutation consisting only of fixed points and
//    pairwise transpositions.
//
//    An involution is its own inverse permutation.
//
//  Recursion:
//
//    S(0) = 1
//    S(1) = 1
//    S(N) = S(N-1) + (N-1) * S(N-2)
//
//  First values:
//
//     N         S(N)
//     0           1
//     1           1
//     2           2
//     3           4
//     4          10
//     5          26
//     6          76
//     7         232
//     8         764
//     9        2620
//    10        9496
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Output, int S[N+1], the number of involutions of 0, 1, 2, ... N
//    objects.
//
{
  int i;

  if ( n < 0 )
  {
    return;
  }

  s[0] = 1;

  if ( n <= 0 )
  {
    return;
  }

  s[1] = 1;

  for ( i = 2; i <= n; i++ )
  {
    s[i] = s[i-1] + ( i - 1 ) * s[i-2];
  }

  return;
}
//****************************************************************************80

void jfrac_to_rfrac ( int m, double r[], double s[], double p[], double q[] )

//****************************************************************************80
//
//  Purpose:
//
//    JFRAC_TO_RFRAC converts a J-fraction into a rational polynomial fraction.
//
//  Discussion:
//
//    The routine accepts a J-fraction:
//
//        R(1) / ( X + S(1)
//      + R(2) / ( X + S(2)
//      + R(3) / ...
//      + R(M) / ( X + S(M) )... ))
//
//    and returns the equivalent rational polynomial fraction:
//
//      P(1) + P(2) * X + ... + P(M) * X**(M-1)
//      -------------------------------------------------------
//      Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
//    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christop Witzgall.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
//    John Rice, Henry Thatcher, Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968.
//
//  Parameters:
//
//    Input, int M, defines the number of P, R, and S
//    coefficients, and is one less than the number of Q
//    coefficients.
//
//    Input, double R[M], S[M], the coefficients defining the J-fraction.
//
//    Output, double P[M], Q[M+1], the coefficients defining the rational
//    polynomial fraction.  The algorithm used normalizes the coefficients
//    so that Q[M+1] = 1.0.
//
{
  double *a;
  double *b;
  int i;
  int k;

  a = new double[m*m];
  b = new double[m*m];

  a[0+0*m] = r[0];
  b[0+0*m] = s[0];

  if ( 1 < m )
  {
    for ( k = 1; k < m; k++ )
    {
      a[k+k*m] = r[0];
      b[k+k*m] = b[k-1+(k-1)*m] + s[k];
    }

    a[0+1*m] = r[0] * s[1];
    b[0+1*m] = r[1] + s[0] * s[1];

    for ( k = 2; k < m; k++ )
    {
      a[0  +k*m] = s[k] * a[0+(k-1)*m]   + r[k] * a[0+(k-2)*m];
      a[k-1+k*m] =        a[k-2+(k-1)*m] + s[k] * r[0];
      b[0  +k*m] = s[k] * b[0+(k-1)*m]   + r[k] * b[0+(k-2)*m];
      b[k-1+k*m] =        b[k-2+(k-1)*m] + s[k] * b[k-1+(k-1)*m] + r[k];
    }

    for ( k = 2; k < m; k++ )
    {
      for ( i = 1; i < k-1; i++ )
      {
        a[i+k*m] = a[i-1+(k-1)*m] + s[k] * a[i+(k-1)*m] 
                                  + r[k] * a[i+(k-2)*m];
        b[i+k*m] = b[i-1+(k-1)*m] + s[k] * b[i+(k-1)*m] 
                                  + r[k] * b[i+(k-2)*m];
      }
    }

  }

  for ( i = 0; i < m; i++ )
  {
    p[i] = a[i+(m-1)*m];
  }

  for ( i = 0; i < m; i++ )
  {
    q[i] = b[i+(m-1)*m];
  }

  q[m] = 1.0;

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

int josephus ( int n, int m, int k )

//****************************************************************************80
//
//  Purpose:
//
//    JOSEPHUS returns the position X of the K-th man to be executed.
//
//  Discussion:
//
//    The classic Josephus problem concerns a circle of 41 men.
//    Every third man is killed and removed from the circle.  Counting
//    and executing continues until all are dead.  Where was the last
//    survivor sitting?
//
//    Note that the first person killed was sitting in the third position.
//    Moreover, when we get down to 2 people, and we need to count the
//    "third" one, we just do the obvious thing, which is to keep counting
//    around the circle until our count is completed.
//
//    The process may be regarded as generating a permutation of
//    the integers from 1 to N.  The permutation would be the execution
//    list, that is, the list of the executed men, by position number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//     John Burkardt.
//
//  Reference:
//
//    W W Rouse Ball,
//    Mathematical Recreations and Essays,
//    Macmillan, 1962, pages 32-36.
//
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 1, Fundamental Algorithms,
//    Addison Wesley, 1968, pages 158-159.
//
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 3, Sorting and Searching,
//    Addison Wesley, 1968, pages 18-19.
//
//  Parameters:
//
//    Input, int N, the number of men.
//    N must be positive.
//
//    Input, int M, the counting index.
//    M must not be zero.  Ordinarily, M is positive, and no greater than N.
//
//    Input, int K, the index of the executed man of interest.
//    K must be between 1 and N.
//
//    Output, int JOSEPHUS, the position of the K-th man.
//    The value will be between 1 and N.
//
{
  int m2;
  int x;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "JOSEPHUS - Fatal error!\n";
    cerr << "  N <= 0.\n";
    exit ( 1 );
  }

  if ( m == 0 )
  {
    cerr << "\n";
    cerr << "JOSEPHUS - Fatal error!\n";
    cerr << "  M = 0.\n";
    exit ( 1 );
  }

  if ( k <= 0 || n < k )
  {
    cerr << "\n";
    cerr << "JOSEPHUS - Fatal error!\n";
    cerr << "  J <= 0 or N < K.\n";
    exit ( 1 );
  }
//
//  In case M is bigger than N, or negative, get the
//  equivalent positive value between 1 and N.
//  You can skip this operation if 1 <= M <= N.
//
  m2 = i4_modp ( m, n );

  x = k * m2;

  while ( n < x )
  {
    x = ( m2 * ( x - n ) - 1 ) / ( m2 - 1 );
  }

  return x;
}
//****************************************************************************80

void ksub_next ( int n, int k, int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_NEXT generates the subsets of size K from a set of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, the desired size of the subsets.  K must
//    be between 0 and N.
//
//    Output, int A[K].  A[I] is the I-th element of the
//    subset.  Thus A[I] will be an integer between 1 and N.
//    Note that the routine will return the values in A
//    in sorted order: 1 <= A[0] < A[1] < ... < A[K-1] <= N
//
//    Input/output, bool &MORE.  Set MORE = FALSE before first call
//    for a new sequence of subsets.  It then is set and remains
//    TRUE as long as the subset computed on this call is not the
//    final one.  When the final subset is computed, MORE is set to
//    FALSE as a signal that the computation is done.
//
{
  int j;
  static int m = 0;
  static int m2 = 0;

  if ( k < 0 || n < k )
  {
    cerr << "\n";
    cerr << "KSUB_NEXT - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  K = " << k << "\n";
    cerr << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  if ( !more )
  {
    m2 = 0;
    m = k;
  }
  else
  {
    if ( m2 < n-m )
    {
      m = 0;
    }
    m = m + 1;
    m2 = a[k-m];
  }

  for ( j = 1; j <= m; j++ )
  {
    a[k+j-m-1] = m2 + j;
  }

  more = ( a[0] != (n-k+1) );

  return;
}
//****************************************************************************80

void ksub_next2 ( int n, int k, int a[], int &in, int &iout )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_NEXT2 generates the subsets of size K from a set of size N.
//
//  Discussion:
//
//    This routine uses the revolving door method.  It has no "memory".
//    It simply calculates the successor of the input set,
//    and will start from the beginning after the last set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//    N must be positive.
//
//    Input, int K, the size of the desired subset.  K must be
//    between 0 and N.
//
//    Input/output, int A[K].  On input, the user must
//    supply a subset of size K in A.  That is, A must
//    contain K unique numbers, in order, between 1 and N.  On
//    output, A(I) is the I-th element of the output subset.
//    The output array is also in sorted order.
//
//    Output, int &IN, the element of the output subset which
//    was not in the input set.  Each new subset differs from the
//    last one by adding one element and deleting another.
//
//    Output, int &IOUT, the element of the input subset which
//    is not in the output subset.
//
{
  int j;
  int m;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "KSUB_NEXT2 - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  but 0 < N is required!\n";
    exit ( 1 );
  }

  if ( k < 0 || n < k )
  {
    cerr << "\n";
    cerr << "KSUB_NEXT2 - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  K = " << k << "\n";
    cerr << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  j = 0;

  for ( ; ; )
  {
    if ( 0 < j || ( k % 2 ) == 0 )
    {
      j = j + 1;

      if ( k < j )
      {
        a[k-1] = k;
        in = k;
        iout = n;
        return;
      }

      if ( a[j-1] != j )
      {
        iout = a[j-1];
        in = iout - 1;
        a[j-1] = in;

        if ( j != 1 )
        {
          in = j - 1;
          a[j-2] = in;
        }

        return;

      }

    }

    j = j + 1;
    m = n;

    if ( j < k )
    {
      m = a[j] - 1;
    }

    if ( m != a[j-1] )
    {
      break;
    }

  }

  in = a[j-1] + 1;
  a[j-1] = in;
  iout = in - 1;

  if ( j != 1 )
  {
    a[j-2] = iout;
    iout = j - 1;
  }

  return;
}
//****************************************************************************80

void ksub_next3 ( int n, int k, int a[], bool &more, int &in, int &iout )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_NEXT3 generates the subsets of size K from a set of size N.
//
//  Discussion:
//
//    The routine uses the revolving door method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//    N must be positive.
//
//    Input, int K, the size of the desired subsets.  K must be
//    between 0 and N.
//
//    Output, int A[K].  A(I) is the I-th element of the
//    output subset.  The elements of A are sorted.
//
//    Input/output, bool &MORE.  On first call, set MORE = FALSE
//    to signal the beginning.  MORE will be set to TRUE, and on
//    each call, the routine will return another K-subset.
//    Finally, when the last subset has been returned,
//    MORE will be set FALSE and you may stop calling.
//
//    Output, int &IN, the element of the output subset which
//    was not in the input set.  Each new subset differs from the
//    last one by adding one element and deleting another.  IN is not
//    defined the first time that the routine returns, and is
//    set to zero.
//
//    Output, int &IOUT, the element of the input subset which is
//    not in the output subset.  IOUT is not defined the first time
//    the routine returns, and is set to zero.
//
{
  int i;
  int j;
  int m;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "KSUB_NEXT3 - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  but 0 < N is required!\n";
    exit ( 1 );
  }

  if ( k < 0 || n < k )
  {
    cerr << "\n";
    cerr << "KSUB_NEXT3 - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  K = " << k << "\n";
    cerr << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  if ( !more )
  {
    in = 0;
    iout = 0;
    i4vec_indicator ( k, a );
    more = ( k != n );
    return;
  }

  j = 0;

  for ( ; ; )
  {
    if ( 0 < j || ( k % 2 ) == 0 )
    {
      j = j + 1;

      if ( a[j-1] != j )
      {
        iout = a[j-1];
        in = iout - 1;
        a[j-1] = in;

        if ( j != 1 )
        {
          in = j - 1;
          a[j-2] = in;
        }

        if ( k != 1 )
        {
          more = ( a[k-2] == k-1 );
        }

        more = ( !more ) || ( a[k-1] != n );

        return;

      }

    }

    j = j + 1;
    m = n;

    if ( j < k )
    {
      m = a[j] - 1;
    }

    if ( m != a[j-1] )
    {
      break;
    }

  }

  in = a[j-1] + 1;
  a[j-1] = in;
  iout = in - 1;

  if ( j != 1 )
  {
    a[j-2] = iout;
    iout = j - 1;
  }

  if ( k != 1 )
  {
    more = ( a[k-2] == k-1 );
  }

  more = ( !more ) || ( a[k-1] != n );

  return;
}
//****************************************************************************80

void ksub_next4 ( int n, int k, int a[], bool &done )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_NEXT4 generates the subsets of size K from a set of size N.
//
//  Discussion:
//
//    The subsets are generated one at a time.
//
//    The routine should be used by setting DONE to TRUE, and then calling
//    repeatedly.  Each call returns with DONE equal to FALSE, the array
//    A contains information defining a new subset.  When DONE returns
//    equal to TRUE, there are no more subsets.
//
//    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such subsets.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the entire set.
//
//    Input, int K, the size of the desired subset.  K must be
//    between 0 and N.
//
//    Input/output, int A[K], contains information about
//    the subsets.  On the first call with DONE = TRUE, the input contents
//    of A don't matter.  Thereafter, the input value of A
//    should be the same as the output value of the previous call.
//    In other words, leave the array alone!
//    On output, as long as DONE is returned FALSE, A contains
//    information defining a subset of K elements of a set of N elements.
//    In other words, A will contain K distinct numbers (in order)
//    between 1 and N.
//
//    Input/output, bool &DONE.
//    On the first call, DONE is an input quantity with a value
//    of TRUE which tells the program to initialize data and
//    return the first subset.
//    On return, DONE is an output quantity that is TRUE as long as
//    the routine is returning another subset, and FALSE when
//    there are no more.
//
{
  int i;
  int j;
  int jsave;

  if ( k < 0 || n < k )
  {
    cerr << "\n";
    cerr << "KSUB_NEXT4 - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  K = " << k << "\n";
    cerr << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }
//
//  First call:
//
  if ( done )
  {
    i4vec_indicator ( k, a );

    if ( 0 < n )
    {
      done = false;
    }
    else
    {
      done = true;
    }
  }
//
//  Next call.
//
  else
  {
    if ( a[0] < n-k+1 )
    {
      done = false;

      jsave = k-1;

      for ( j = 0; j < k-1; j++ )
      {
        if ( a[j] + 1 < a[j+1] )
        {
          jsave = j;
          break;
        }

      }

      i4vec_indicator ( jsave, a );
      a[jsave] = a[jsave] + 1;
    }
    else
    {
      done = true;
    }
  }
  return;
}
//****************************************************************************80

void ksub_random ( int n, int k, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_RANDOM selects a random subset of size K from a set of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, number of elements in desired subsets.  K must
//    be between 0 and N.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int A[K].  A(I) is the I-th element of the
//    output set.  The elements of A are in order.
//
{
  int i;
  int ids;
  int ihi;
  int ip;
  int ir;
  int is;
  int ix;
  int l;
  int ll;
  int m;
  int m0;
  double r;

  if ( k < 0 )
  {
    cerr << "\n";
    cerr << "KSUB_RANDOM - Fatal error!\n";
    cerr << "  K = " << k << "\n";
    cerr << "  but 0 <= K is required!\n";
    exit ( 1 );
  }
  else if ( n < k )
  {
    cerr << "\n";
    cerr << "KSUB_RANDOM - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  K = " << k << "\n";
    cerr << "  K <= N is required!\n";
    exit ( 1 );
  }

  if ( k == 0 )
  {
    return;
  }

  for ( i = 1; i <= k; i++ )
  {
    a[i-1] = ( ( i - 1 ) * n ) / k;
  }

  for ( i = 1; i <= k; i++ )
  {
    for ( ; ; )
    {
      ix = i4_uniform ( 1, n, seed );

      l = 1 + ( ix * k - 1 ) / n;

      if ( a[l-1] < ix )
      {
        break;
      }

    }

    a[l-1] = a[l-1] + 1;

  }

  ip = 0;
  is = k;

  for ( i = 1; i <= k; i++ )
  {
    m = a[i-1];
    a[i-1] = 0;

    if ( m != ( ( i - 1 ) * n ) / k )
    {
      ip = ip + 1;
      a[ip-1] = m;
    }

  }

  ihi = ip;

  for ( i = 1; i <= ihi; i++ )
  {
    ip = ihi + 1 - i;
    l = 1 + ( a[ip-1] * k - 1 ) / n;
    ids = a[ip-1] - ( ( l - 1 ) * n ) / k;
    a[ip-1] = 0;
    a[is-1] = l;
    is = is - ids;
  }

  for ( ll = 1; ll <= k; ll++ )
  {
    l = k + 1 - ll;

    if ( a[l-1] != 0 )
    {
      ir = l;
      m0 = 1 + ( ( a[l-1] - 1 ) * n ) / k;
      m = ( a[l-1] * n ) / k - m0 + 1;
    }

    ix = i4_uniform ( m0, m0 + m - 1, seed );

    i = l + 1;

    while ( i <= ir )
    {
      if ( ix < a[i-1] )
      {
        break;
      }

      ix = ix + 1;
      a[i-2] = a[i-1];
      i = i + 1;
    }
    a[i-2] = ix;
    m = m - 1;
  }
  return;
}
//****************************************************************************80

void ksub_random2 ( int n, int k, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_RANDOM2 selects a random subset of size K from a set of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, number of elements in desired subsets.  K must
//    be between 0 and N.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int A[K].  A(I) is the I-th element of the
//    output set.  The elements of A are in order.
//
{
  int available;
  int candidate;
  int have;
  int need;
  double r;

  if ( k < 0 || n < k )
  {
    cerr << "\n";
    cerr << "KSUB_RANDOM2 - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  K = " << k << "\n";
    cerr << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  if ( k == 0 )
  {
    return;
  }

  need = k;
  have = 0;
  available = n;
  candidate = 0;

  for ( ; ; )
  {
    candidate = candidate + 1;

    r = r8_uniform_01 ( seed );

    if ( r * ( double ) available <= ( double ) need )
    {
      need = need - 1;
      a[have] = candidate;
      have = have + 1;

      if ( need <= 0 )
      {
        break;
      }
    }
    available = available - 1;
  }
  return;
}
//****************************************************************************80

void ksub_random3 ( int n, int k, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_RANDOM3 selects a random subset of size K from a set of size N.
//
//  Discussion:
//
//    This routine uses Floyd's algorithm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, number of elements in desired subsets.  K must
//    be between 0 and N.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int A[N].  I is an element of the subset
//    if A(I) = 1, and I is not an element if A(I)=0.
//
{
  int i;
  int j;
  double r;

  if ( k < 0 || n < k )
  {
    cerr << "\n";
    cerr << "KSUB_RANDOM3 - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  K = " << k << "\n";
    cerr << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }

  if ( k == 0 )
  {
    return;
  }

  for ( i = n - k + 1; i <= n; i++ )
  {
    j = i4_uniform ( 1, i, seed );

    if ( a[j-1] == 0 )
    {
      a[j-1] = 1;
    }
    else
    {
      a[i-1] = 1;
    }
  }
  return;
}
//****************************************************************************80

void ksub_random4 ( int n, int k, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_RANDOM4 selects a random subset of size K from a set of size N.
//
//  Discussion:
//
//    This routine is somewhat impractical for the given problem, but
//    it is included for comparison, because it is an interesting
//    approach that is superior for certain applications.
//
//    The approach is mainly interesting because it is "incremental";
//    it proceeds by considering every element of the set, and does not
//    need to know how many elements there are.
//
//    This makes this approach ideal for certain cases, such as the
//    need to pick 5 lines at random from a text file of unknown length,
//    or to choose 6 people who call a certain telephone number on a
//    given day.  Using this technique, it is possible to make the
//    selection so that, whenever the input stops, a valid uniformly
//    random subset has been chosen.
//
//    Obviously, if the number of items is known in advance, and
//    it is easy to extract K items directly, there is no need for
//    this approach, and it is less efficient since, among other costs,
//    it has to generate a random number for each item, and make an
//    acceptance/rejection test.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Tom Christiansen, Nathan Torkington,
//    "8.6: Picking a Random Line from a File",
//    Perl Cookbook, pages 284-285,
//    O'Reilly, 1999.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, number of elements in desired subsets.  K must
//    be between 0 and N.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int A[K], contains the indices of the selected items.
//
{
  int i;
  int next;
  double r;

  next = 0;
//
//  Here, we use a WHILE to suggest that the algorithm
//  proceeds to the next item, without knowing how many items
//  there are in total.
//
//  Note that this is really the only place where N occurs,
//  so other termination criteria could be used, and we really
//  don't need to know the value of N!
//
  while ( next < n )
  {
    next = next + 1;

    if ( next <= k )
    {
      i = next;
      a[i-1] = next;
    }
    else
    {
      r = r8_uniform_01 ( seed );

      if ( r * ( double ) next <= ( double ) k )
      {
        i = i4_uniform ( 1, k, seed );
        a[i-1] = next;
      }
    }
  }
  return;
}
//****************************************************************************80

int *ksub_random5 ( int n, int k, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_RANDOM5 selects a random subset of size K from a set of size N.
//
//  Discussion:
//
//    Consider the set A = 1, 2, 3, ... N.  
//    Choose a random index I1 between 1 and N, and swap items A(1) and A(I1).
//    Choose a random index I2 between 2 and N, and swap items A(2) and A(I2).
//    repeat K times.
//    A(1:K) is your random K-subset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets
//    are drawn.
//
//    Input, int K, number of elements in desired subsets.
//    1 <= K <= N.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, int KSUB_RANDOM5[K], the indices of the randomly
//    chosen elements.  These are 1-based indices.
//
{
  int *a;
  int *b;
  static int base = 1;
  int i;
  int j;
  int t;
//
//  Let B index the set.
//
  b = new int[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = i + base;
  }
//
//  Choose item 1 from N things,
//  choose item 2 from N-1 things,
//  choose item K from N-K+1 things.
//
  for ( i = 0; i < k; i++)
  {
    j = i4_uniform ( i, n - 1, seed );
    t    = b[i];
    b[i] = b[j];
    b[j] = t;
  }
//
//  Copy the first K elements.
//
  a = new int[k];

  for ( i = 0; i < k; i++ )
  {
    a[i] = b[i];
  }
  delete [] b;
//
//  Put the elements in ascending order.
//
  i4vec_sort_heap_a ( k, a );

  return a;
}
//****************************************************************************80

void ksub_rank ( int k, int a[], int &rank )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_RANK computes the rank of a subset of an N set.
//
//  Discussion:
//
//    The routine accepts an array representing a subset of size K from a set
//    of size N, and returns the rank (or order) of that subset.  
//
//    This is the same order in which routine KSUB_NEXT2 would produce that subset.
//
//    Note the value of N is not input, and is not, in fact,
//    needed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int K, the number of elements in the subset.
//
//    Input, int A[K], contains K distinct numbers between
//    1 and N, in order.
//
//    Output, int &RANK, the rank of this subset.
//
{
  int i;
  int iprod;
  int j;

  rank = 0;

  for ( i = 1; i <= k; i++ )
  {
    iprod = 1;

    for ( j = i+1; j <= a[i-1]-1; j++ )
    {
      iprod = iprod * j;
    }

    for ( j = 1; j <= a[i-1]-i-1; j++ )
    {
      iprod = iprod / j;
    }

    if ( a[i-1] == 1 )
    {
      iprod = 0;
    }
    rank = rank + iprod;
  }

  rank = rank + 1;

  return; 
}
//****************************************************************************80

void ksub_unrank ( int k, int rank, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_UNRANK returns the subset of a given rank.
//
//  Discussion:
//
//    The routine is given a rank and returns the corresponding subset of K
//    elements of a set of N elements.  
//
//    It uses the same ranking that KSUB_NEXT2 uses to generate all the subsets 
//    one at a time.  
//
//    Note that the value of N itself is not input, nor is it needed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int K, the number of elements in the subset.
//
//    Input, int RANK, the rank of the desired subset.
//    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such
//    subsets, so RANK must be between 1 and that value.
//
//    Output, int A[K], K distinct integers in order between
//    1 and N, which define the subset.
//
{
  int i;
  int ii;
  int ip;
  int iprod;
  int jrank;

  jrank = rank - 1;

  for ( i = k; 1 <= i; i-- )
  {
    ip = i - 1;
    iprod = 1;

    for ( ; ; )
    {
      ip = ip + 1;

      if ( ip != i )
      {
        iprod = ( ip * iprod ) / ( ip - i );
      }

      if ( jrank < iprod )
      {
        break;
      }
    }

    if ( ip != i )
    {
      iprod = ( ( ip - i ) * iprod ) / ip;
    }

    jrank = jrank - iprod;
    a[i-1] = ip;
  }

  return;
}
//****************************************************************************80

void lvec_next ( int n, bool lvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    LVEC_NEXT generates the next logical vector.
//
//  Discussion:
//
//    In the following discussion, we will let '0' stand for FALSE and
//    '1' for TRUE.
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
//    31 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input/output, bool LVEC[N], on output, the successor to the
//    input vector.  
//
{
  int i;

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( !lvec[i] )
    {
      lvec[i] = true;
      return;
    }
    lvec[i] = false;
  }
  return;
}
//****************************************************************************80

void matrix_product_opt ( int n, int rank[], int &cost, int order[] )

//****************************************************************************80
//
//  Purpose:
//
//    MATRIX_PRODUCT_OPT determines the optimal cost of a matrix product.
//
//  Discussion:
//
//    The cost of multiplying an LxM matrix by an M by N matrix is
//    assessed as L*M*N.
//
//    Any particular order of multiplying a set of N matrices is equivalent
//    to parenthesizing an expression of N objects.
//
//    The actual number of ways of parenthesizing an expression
//    of N objects is C(N), the N-th Catalan number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Sedgewick,
//    Algorithms,
//    Addison-Wesley, 1984, pages 486-489.
//
//  Parameters:
//
//    Input, int N, the number of matrices to be multiplied.
//
//    Input, int RANK[N+1], the rank information for the matrices.
//    Matrix I has RANK[I] rows and RANK[I+1] columns.
//
//    Output, int &COST, the cost of the multiplication if the optimal
//    order is used.
//
//    Output, int ORDER[N-1], indicates the order in which the N-1
//    multiplications are to be carried out.  ORDER[0] is the first
//    multiplication to do, and so on.
//
{
# define STACK_MAX 100

  int *best;
  int *cost2;
  int cost3;
  int i;
  int i_inf;
  int i1;
  int i2;
  int i3;
  int j;
  int k;
  int stack[STACK_MAX];
  int stack_num;
  int step;
//
//  Initialize the cost matrix.
//
  best = new int[n*n];
  cost2 = new int[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j <= i; j++ )
    {
      cost2[i+j*n] = 0;
    }
    for ( j = i+1; j < n; j++ )
    {
      cost2[i+j*n] = i4_huge ( );
    }
  }
//
//  Initialize the BEST matrix.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      best[i+j*n] = 0;
    }
  }
//
//  Compute the cost and best matrices.
//
  for ( j = 1; j <= n-1; j++ )
  {
    for ( i = 1; i <= n-j; i++ )
    {
      for ( k = i+1; k <= i+j; k++ )
      {
        cost3 = cost2[i-1+(k-2)*n] + cost2[k-1+(i+j-1)*n] 
          + rank[i-1] * rank[k-1] * rank[i+j];

        if ( cost3 < cost2[i-1+(i+j-1)*n] )
        {
          cost2[i-1+(i+j-1)*n] = cost3;
          best[i-1+(i+j-1)*n] = k;
        }
      }
    }
  }
//
//  Pick off the optimal cost.
//
  cost = cost2[0+(n-1)*n];
//
//  Backtrack to determine the optimal order.
//
  stack_num = 0;

  i1 = 1;
  i2 = n;

  if ( i1+1 < i2 )
  {
    stack[stack_num] = i1;
    stack_num = stack_num + 1;
    stack[stack_num] = i2;
    stack_num = stack_num + 1;
  }

  step = n - 1;
//
//  Take an item off the stack.
//
  while ( 0 < stack_num )
  {
    stack_num = stack_num - 1;
    i3 = stack[stack_num];
    stack_num = stack_num - 1;
    i1 = stack[stack_num];

    i2 = best[i1-1+(i3-1)*n];

    step = step - 1;
    order[step] = i2 - 1;
//
//  The left chunk is matrices (I1...I2-1)
//
    if ( i1 == i2-1 )
    {
    }
    else if ( i1+1 == i2-1 )
    {
      step = step - 1;
      order[step] = i2 - 2;
    }
    else
    {
      stack[stack_num] = i1;
      stack_num = stack_num + 1;
      stack[stack_num] = i2 - 1;
      stack_num = stack_num + 1;
    }
//
//  The right chunk is matrices (I2...I3)
//
    if ( i2 == i3 )
    {
    }
    else if ( i2+1 == i3 )
    {
      step = step - 1;
      order[step] = i2;
    }
    else
    {
      stack[stack_num] = i2;
      stack_num = stack_num + 1;
      stack[stack_num] = i3;
      stack_num = stack_num + 1;
    }

  }
  delete [] best;
  delete [] cost2;

  return;
# undef STACK_MAX
}
//****************************************************************************80

void moebius_matrix ( int n, int a[], int mu[] )

//****************************************************************************80
//
//  Purpose:
//
//    MOEBIUS_MATRIX finds the Moebius matrix from a covering relation.
//
//  Discussion:
//
//    This routine can be called with A and MU being the same matrix.
//    The routine will correctly compute the Moebius matrix, which
//    will, in this case, overwrite the input matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, number of elements in the partially ordered set.
//
//    Input, int A[N*N].  A(I,J) = 1 if I is covered by J,
//    0 otherwise.
//
//    Output, int MU[N*N], the Moebius matrix as computed by the routine.
//
{
  int i;
  int j;
  int *p;

  p = new int[n];
//
//  Compute a reordering of the elements of the partially ordered matrix.
//
  triang ( n, a, p );
//
//  Copy the matrix.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      mu[i+j*n] = a[i+j*n];
    }
  }
//
//  Apply the reordering to MU.
//
  i4mat_perm2 ( n, n, mu, p, p );
//
//  Negate the (strict) upper triangular elements of MU.
//
  for ( i = 0; i < n-1; i++ )
  {
    for ( j = i+1; j < n; j++ )
    {
      mu[i+j*n] = - mu[i+j*n];
    }
  }
//
//  Compute the inverse of MU.
//
  i4mat_u1_inverse ( n, mu, mu );
//
//  All nonzero elements are reset to 1.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = i; j < n; j++ )
    {
      if ( mu[i+j*n] != 0 )
      {
        mu[i+j*n] = 1;
      }
    }
  }
//
//  Invert the matrix again.
//
  i4mat_u1_inverse ( n, mu, mu );
//
//  Compute the inverse permutation.
//
  perm_inverse ( n, p );
//
//  Unpermute the rows and columns of MU.
//
  i4mat_perm2 ( n, n, mu, p, p );

  delete [] p;

  return;
}
//****************************************************************************80

int monomial_count ( int degree_max, int dim )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_COUNT counts the number of monomials up to a given degree.
//
//  Discussion:
//
//    In 3D, there are 10 monomials of degree 3 or less:
//
//    Degree  Count  List
//    ------  -----  ----
//         0      1  1
//         1      3  x y z
//         2      6  xx xy xz yy yz zz
//         3     10  xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
//
//    Total      20
//
//    The formula is 
//
//      COUNTS(DEGREE,DIM) = (DIM-1+DEGREE)! / (DIM-1)! / DEGREE!
//
//      TOTAL              = (DIM  +DEGREE)! / (DIM)!   / DEGREE!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum degree.
//
//    Input, int DIM, the spatial dimension.
//
//    Output, int MONOMIAL_COUNT, the total number of monomials
//    of degrees 0 through DEGREE_MAX.
//
{
  int bot;
  int top;
  int total;

  total = 1;

  if ( degree_max < dim )
  {
    top = dim + 1;
    for ( bot = 1; bot <= degree_max; bot++ )
    {
      total = ( total * top ) / bot;
      top = top + 1;
    }
  }
  else
  {
    top = degree_max + 1;
    for ( bot = 1; bot <= dim; bot++ )
    {
      total = ( total * top ) / bot;
      top = top + 1;
    }
  }
  return total;
}
//****************************************************************************80

int *monomial_counts ( int degree_max, int dim )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_COUNTS counts the number of monomials up to a given degree.
//
//  Discussion:
//
//    In 3D, there are 10 monomials of degree 3 or less:
//
//    Degree  Count  List
//    ------  -----  ----
//         0      1  1
//         1      3  x y z
//         2      6  xx xy xz yy yz zz
//         3     10  xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
//
//    Total      20
//
//    The formula is 
//
//      COUNTS(DEGREE,DIM) = (DIM-1+DEGREE)! / (DIM-1)! / DEGREE!
//
//      TOTAL              = (DIM  +DEGREE)! / (DIM)!   / DEGREE!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum degree.
//
//    Input, int DIM, the spatial dimension.
//
//    Output, int MONOMIAL_COUNTS[DEGREE_MAX+1], the number of
//    monomials of each degree.
//
{
  int *counts;
  int degree;

  counts = new int[degree_max+1];

  degree = 0;
  counts[degree] = 1;

  for ( degree = 1; degree <= degree_max; degree++ )
  {
    counts[degree] = ( counts[degree-1] * ( dim - 1 + degree ) ) / degree;
  }

  return counts;
}
//****************************************************************************80

int morse_thue ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    MORSE_THUE generates a Morse_Thue number.
//
//  Discussion:
//
//    The Morse_Thue sequence can be defined in a number of ways.
//
//    A) Start with the string containing the single letter '0'; then
//       repeatedly apply the replacement rules '0' -> '01' and
//       '1' -> '10' to the letters of the string.  The Morse_Thue sequence
//       is the resulting letter sequence.
//
//    B) Starting with the string containing the single letter '0',
//       repeatedly append the binary complement of the string to itself.
//       Thus, '0' becomes '0' + '1' or '01', then '01' becomes
//       '01' + '10', which becomes '0110' + '1001', and so on.
//
//    C) Starting with I = 0, the I-th Morse-Thue number is determined
//       by taking the binary representation of I, adding the digits,
//       and computing the remainder modulo 2.
//
//  Example:
//
//     I  binary   S
//    --  ------  --
//     0       0   0
//     1       1   1
//     2      10   1
//     3      11   0
//     4     100   1
//     5     101   0
//     6     110   0
//     7     111   1
//     8    1000   1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the index of the Morse-Thue number.
//    Normally, I is 0 or greater, but any value is allowed.
//
//    Output, int MORSE_THUE, the Morse-Thue number of index I.
//
{
# define NBITS 32

  int b[NBITS];
  int s;

  i = abs ( i );
//
//  Expand I into binary form.
//
  ui4_to_ubvec ( i, NBITS, b );
//
//  Sum the 1's in the binary representation.
//
  s = 0;
  for ( i = 0; i < NBITS; i++ )
  {
    s = s + b[i];
  }
//
//  Take the value modulo 2.
//
  s = s % 2;

  return s;
# undef NBITS
}
//****************************************************************************80

int multinomial_coef1 ( int nfactor, int factor[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_COEF1 computes a multinomial coefficient.
//
//  Discussion:
//
//    The multinomial coefficient is a generalization of the binomial
//    coefficient.  It may be interpreted as the number of combinations of
//    N objects, where FACTOR(1) objects are indistinguishable of type 1,
//    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
//    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
//
//    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
//
//    The log of the gamma function is used, to avoid overflow.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NFACTOR, the number of factors.
//
//    Input, int FACTOR[NFACTOR], contains the factors.
//    0 <= FACTOR(I)
//
//    Output, int MULTINOMIAL_COEF1, the value of the multinomial coefficient.
//
{
  double arg;
  double fack;
  double facn;
  int i;
  int n;
  int value;
//
//  Each factor must be nonnegative.
//
  for ( i = 0; i < nfactor; i++ )
  {
    if ( factor[i] < 0 )
    {
      cerr << "\n";
      cerr << "MULTINOMIAL_COEF1 - Fatal error\n";
      cerr << "  Factor " << i << " = " << factor[i] << "\n";
      cerr << "  But this value must be nonnegative.\n";
      exit ( 1 );
    }
  }
//
//  The factors sum to N.
//
  n = i4vec_sum ( nfactor, factor );

  arg = ( double ) ( n + 1 );
  facn = r8_gamma_log ( arg );

  for ( i = 0; i < nfactor; i++ )
  {
    arg = ( double ) ( factor[i] + 1 );
    fack = r8_gamma_log ( arg );
    facn = facn - fack;
  }

  value = r8_nint ( exp ( facn ) );

  return value;
}
//****************************************************************************80

int multinomial_coef2 ( int nfactor, int factor[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_COEF2 computes a multinomial coefficient.
//
//  Discussion:
//
//    The multinomial coefficient is a generalization of the binomial
//    coefficient.  It may be interpreted as the number of combinations of
//    N objects, where FACTOR(1) objects are indistinguishable of type 1,
//    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
//    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
//
//    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
//
//    A direct method is used, which should be exact.  However, there
//    is a possibility of intermediate overflow of the result.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NFACTOR, the number of factors.
//
//    Input, int FACTOR[NFACTOR], contains the factors.
//    0 <= FACTOR(I)
//
//    Output, int MULTINOMIAL_COEF2, the value of the multinomial coefficient.
//
{
  int i;
  int j;
  int k;
  int value;
//
//  Each factor must be nonnegative.
//
  for ( i = 0; i < nfactor; i++ )
  {
    if ( factor[i] < 0 )
    {
      cerr << "\n";
      cerr << "MULTINOMIAL_COEF2 - Fatal error!\n";
      cerr << "  Factor " << i << " = " << factor[i] << "\n";
      cerr << "  But this value must be nonnegative.\n";
      exit ( 1 );
    }
  }

  value = 1;
  k = 0;

  for ( i = 0; i < nfactor; i++ )
  {
    for ( j = 1; j <= factor[i]; j++ )
    {
      k = k + 1;
      value = ( value * k ) / j;
    }
  }

  return value;
}
//****************************************************************************80

int multiperm_enum ( int n, int k, int counts[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIPERM_ENUM enumerates multipermutations.
//
//  Discussion:
//
//    A multipermutation is a permutation of objects, some of which are
//    identical.
//
//    While there are 6 permutations of the distinct objects A,B,C, there
//    are only 3 multipermutations of the objects A,B,B.
//
//    In general, there are N! permutations of N distinct objects, but
//    there are N! / ( ( M1! ) ( M2! ) ... ( MK! ) ) multipermutations
//    of N objects, in the case where the N objects consist of K
//    types, with M1 examples of type 1, M2 examples of type 2 and so on,
//    and for which objects of the same type are indistinguishable.
//
//  Example:
//
//    Input:
//
//      N = 5, K = 3, COUNTS = (/ 1, 2, 2 /)
//
//    Output:
//
//      Number = 30
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items in the multipermutation.
//
//    Input, int K, the number of types of items.
//    1 <= K.  Ordinarily, K <= N, but we allow any positive K, because
//    we also allow entries in COUNTS to be 0.
//
//    Input, int COUNTS[K], the number of items of each type.
//    0 <= COUNTS(1:K) <= N and sum ( COUNTS(1:K) ) = N.
//
//    Output, int MULTIPERM_ENUM, the number of multipermutations.
//
{
  int i;
  int j;
  int number;
  int sum;
  int top;

  if ( n < 0 )
  {
    number = -1;
    return number;
  }

  if ( n == 0 )
  {
    number = 1;
    return number;
  }

  if ( k < 1 )
  {
    number = -1;
    return number;
  }

  for ( i = 0; i < k; i++ )
  {
    if ( counts[i] < 0 )
    {
      number = -1;
      return number;
    }
  }

  sum = 0;
  for ( i = 0; i < k; i++ )
  {
    sum = sum + counts[i];
  }
  if ( sum != n )
  {
    number = -1;
    return number;
  }
//
//  Ready for computation.
//  By design, the integer division should never have a remainder.
//
  top = 0;
  number = 1;

  for ( i = 0; i < k; i++ )
  {
    for ( j = 1; j <= counts[i]; j++ )
    {
      top = top + 1;
      number = ( number * top ) / j;
    }
  }

  return number;
}
//****************************************************************************80

void multiperm_next ( int n, int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIPERM_NEXT returns the next multipermutation.
//
//  Discussion:
//
//    To begin the computation, the user must set up the first multipermutation.
//    To compute ALL possible multipermutations, this first permutation should
//    list the values in ascending order.
//
//    The routine will compute, one by one, the next multipermutation,
//    in lexicographical order.  On the call after computing the last 
//    multipermutation, the routine will return MORE = FALSE (and will
//    reset the multipermutation to the FIRST one again.)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items in the multipermutation.
//
//    Input/output, int A[N]; on input, the current multipermutation.
//    On output, the next multipermutation.
//
//    Output, bool &MORE, is TRUE if the next multipermutation
//    was computed, or FALSE if no further multipermutations could
//    be computed.
//
{
  int i;
  int m;
  int temp;
//
//  Step 1:
//  Find M, the last location in A for which A(M) < A(M+1).
//
  m = 0;
  for ( i = 1; i <= n-1; i++ )
  {
    if ( a[i-1] < a[i] )
    {
      m = i;
    }
  }
//
//  Step 2:
//  If no M was found, we've run out of multipermutations.
//
  if ( m == 0 )
  {
    more = false;
    i4vec_sort_heap_a ( n, a );
    return;
  }
  else
  {
    more = true;
  }
//
//  Step 3:
//  Ascending sort A(M+1:N).
//
  if ( m + 1 < n )
  {
    i4vec_sort_heap_a ( n - m, a + m );
  }
//
//  Step 4:
//  Locate the first larger value after A(M).
//
  i = 1;
  for ( ; ; )
  {
    if ( a[m-1] < a[m+i-1] )
    {
      break;
    }
    i = i + 1;
  }
//
//  Step 5:
//  Interchange A(M) and the next larger value.
//
  temp = a[m-1];
  a[m-1] = a[m+i-1];
  a[m+i-1] = temp;

  return;
}
//****************************************************************************80

void network_flow_max ( int nnode, int nedge, int iendpt[], int icpflo[], 
  int source, int sink, int cut[], int node_flow[] )

//****************************************************************************80
//
//  Purpose:
//
//    NETWORK_FLOW_MAX finds the maximal flow and a minimal cut in a network.
//
//  Discussion:
//
//    Apparently, I didn't get around to converting this routine!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2000
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int NNODE, the number of nodes.
//
//    Input, int NEDGE, the number of edges.
//
//    Input/output, int IENDPT[2*NEDGE], the edges of the network,
//    defined as pairs of nodes.  Each edge should be listed TWICE,
//    the second time in reverse order.  On output, the edges have
//    been reordered, and so the columns of IENDPT have been rearranged.
//
//    Input/output, int ICPFLO[2*NEDGE].  Capacities and flows.
//    On input, ICPFLO(1,I) is the capacity of edge I.  On output,
//    ICPFLO(2,I) is the flow on edge I and ICPFLO(1,I) has
//    been rearranged to match the reordering of IENDPT.
//
//    Input, int SOURCE, the designated source node.
//
//    Input, int SINK, the designated sink node.
//
//    Output, int CUT[NNODE].  CUT(I) = 1 if node I is in the
//    minimal cut set, otherwise 0.
//
//    Output, int NODE_FLOW[NNODE].  NODE_FLOW(I) is the value of the flow
//    through node I.
//
{
  return;
}
//****************************************************************************80

unsigned int nim_sum ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    NIM_SUM computes the Nim sum of two integers.
//
//  Discussion:
//
//    If K is the Nim sum of I and J, then each bit of K is the exclusive
//    OR of the corresponding bits of I and J.
//
//  Example:
//
//     I     J     K       I_2        J_2         K_2
//   ----  ----  ----  ----------  ----------  ----------
//      0     0     0           0           0           0
//      1     0     1           1           0           1
//      1     1     0           1           1           0
//      2     7     5          10         111         101
//     11    28    23        1011       11100       10111
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the integers to be Nim-summed.
//
//    Output, unsigned int NIM_SUM, the Nim sum of I and J.
//
{
# define NBITS 32

  int bvec1[NBITS];
  int bvec2[NBITS];
  int bvec3[NBITS];
  unsigned int value;

  ui4_to_ubvec ( i, NBITS, bvec1 );

  ui4_to_ubvec ( j, NBITS, bvec2 );

  bvec_xor ( NBITS, bvec1, bvec2, bvec3 );

  value = ubvec_to_ui4 ( NBITS, bvec3 );

  return value;

# undef NBITS
}
//****************************************************************************80

void padovan ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PADOVAN returns the first N values of the Padovan sequence.
//
//  Discussion:
//
//    The Padovan sequence has the initial values:
//
//      P(0) = 1
//      P(1) = 1
//      P(2) = 1
//
//    and subsequent entries are generated by the recurrence
//
//      P(I+1) = P(I-1) + P(I-2)
//
//  Example:
//
//    0   1
//    1   1
//    2   1
//    3   2
//    4   2
//    5   3
//    6   4
//    7   5
//    8   7
//    9   9
//   10  12
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Stewart,
//    "A Neglected Number",
//    Scientific American, Volume 274, pages 102-102, June 1996.
//
//    Ian Stewart,
//    Math Hysteria,
//    Oxford, 2004.
//
//  Parameters:
//
//    Input, int N, the number of terms.
//
//    Output, int P[N], terms 0 though N-1 of the sequence.
//
{
  int i;

  if ( n < 1 )
  {
    return;
  }

  p[0] = 1;

  if ( n < 2 )
  {
    return;
  }

  p[1] = 1;

  if ( n < 3 )
  {
    return;
  }
 
  p[2] = 1;

  for ( i = 4; i <= n; i++ )
  {
    p[i-1] = p[i-3] + p[i-4];
  }

  return;
}
//****************************************************************************80

void pell_basic ( int d, int &x0, int &y0 )

//****************************************************************************80
//
//  Purpose:
//
//    PELL_BASIC returns the fundamental solution for Pell's basic equation.
//
//  Discussion:
//
//    Pell's equation has the form:
//
//      X**2 - D * Y**2 = 1
//
//    where D is a given non-square integer, and X and Y may be assumed
//    to be positive integers.
//
//  Example:
//
//     D   X0   Y0
//
//     2    3    2
//     3    2    1
//     5    9    4
//     6    5    2
//     7    8    3
//     8    3    1
//    10   19    6
//    11   10    3
//    12    7    2
//    13  649  180
//    14   15    4
//    15    4    1
//    17   33    8
//    18   17    4
//    19  170   39
//    20    9    2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//   John Burkardt
//
//  Reference:
//
//    Mark Herkommer,
//    Number Theory, A Programmer's Guide,
//    McGraw Hill, 1999, pages 294-307
//
//  Parameters:
//
//    Input, int D, the coefficient in Pell's equation.  D should be
//    positive, and not a perfect square.
//
//    Output, int &X0, &Y0, the fundamental or 0'th solution.
//    If X0 = Y0 = 0, then the calculation was canceled because of an error.
//    Both X0 and Y0 will be nonnegative.
//
{
# define TERM_MAX 100

  int b[TERM_MAX+1];
  int i;
  int p;
  int pm1;
  int pm2;
  int q;
  int qm1;
  int qm2;
  int r;
  int term_num;
//
//  Check D.
//
  if ( d <= 0 )
  {
    cerr << "\n";
    cerr << "PELL_BASIC - Fatal error!\n";
    cerr << "  Pell coefficient D <= 0.\n";
    x0 = 0;
    y0 = 0;
    exit ( 1 );
  }

  i4_sqrt ( d, q, r );

  if ( r == 0 )
  {
    cerr << "\n";
    cerr << "PELL_BASIC - Fatal error!\n";
    cerr << "  Pell coefficient is a perfect square.\n";
    x0 = 0;
    y0 = 0;
    exit ( 1 );
  }
//
//  Find the continued fraction representation of sqrt ( D ).
//
  i4_sqrt_cf ( d, TERM_MAX, term_num, b );
//
//  If necessary, go for two periods.
//
  if ( ( term_num % 2 ) == 1 )
  {
    for ( i = term_num+1; i <= 2 * term_num; i++ )
    {
      b[i] = b[i-term_num];
    }
    term_num = 2 * term_num;
  }
//
//  Evaluate the continued fraction using the forward recursion algorithm.
//
  pm2 = 0;
  pm1 = 1;
  qm2 = 1;
  qm1 = 0;

  for ( i = 0; i < term_num; i++ )
  {
    p = b[i] * pm1 + pm2;
    q = b[i] * qm1 + qm2;
    pm2 = pm1;
    pm1 = p;
    qm2 = qm1;
    qm1 = q;
  }
//
//  Get the fundamental solution.
//
  x0 = p;
  y0 = q;

  return;
# undef TERM_MAX
}
//****************************************************************************80

void pell_next ( int d, int x0, int y0, int xn, int yn, int &xnp1, int &ynp1 )

//****************************************************************************80
//
//  Purpose:
//
//    PELL_NEXT returns the next solution of Pell's equation.
//
//  Discussion:
//
//    Pell's equation has the form:
//
//      X**2 - D * Y**2 = 1
//
//    where D is a given non-square integer, and X and Y may be assumed
//    to be positive integers.
//
//    To compute X0, Y0, call PELL_BASIC.
//    To compute X1, Y1, call this routine, with XN and YN set to X0 and Y0.
//    To compute further solutions, call again with X0, Y0 and the previous
//    solution.
//
//  Example:
//
//    ------INPUT--------  --OUTPUT--
//
//    D  X0  Y0   XN   YN  XNP1  YNP1
//
//    2   3   2    3    2    17    12
//    2   3   2   17   12    99    70
//    2   3   2   99   70   577   408
//    2   3   2  577  408  3363  2378
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//   John Burkardt
//
//  Reference:
//
//    Mark Herkommer,
//    Number Theory, A Programmer's Guide,
//    McGraw Hill, 1999, pages 294-307
//
//  Parameters:
//
//    Input, int D, the coefficient in Pell's equation.
//
//    Input, int X0, Y0, the fundamental or 0'th solution.
//
//    Input, int XN, YN, the N-th solution.
//
//    Output, int &XNP1, &YNP1, the N+1-th solution.
//
{
  xnp1 = x0 * xn + d * y0 * yn;
  ynp1 = x0 * yn +     y0 * xn;

  return;
}
//****************************************************************************80

int pent_enum ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PENT_ENUM computes the N-th pentagonal number.
//
//  Discussion:
//
//    The pentagonal number P(N) counts the number of dots in a figure of
//    N nested pentagons.  The pentagonal numbers are defined for both
//    positive and negative N.
//
//  First values:
//
//    N   P
//
//   -5   40
//   -4   26
//   -3   15
//   -2    7
//   -1    2
//    0    0
//    1    1
//    2    5
//    3   12
//    4   22
//    5   35
//    6   51
//    7   70
//    8   92
//    9  117
//   10  145
//
//    P(N) = ( N * ( 3 * N - 1 ) ) / 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the index of the pentagonal number desired.
//
//    Output, int PENT_ENUM, the value of the N-th pentagonal number.
//
{
  int p;

  p = ( n * ( 3 * n - 1 ) ) / 2;

  return p;
}
//****************************************************************************80

void perm_ascend ( int n, int a[], int &length, int sub[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_ASCEND computes the longest ascending subsequence of permutation.
//
//  Discussion:
//
//    Although this routine is intended to be applied to a permutation,
//    it will work just as well for an arbitrary vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the permutation.
//
//    Input, int A[N], the permutation to be examined.
//
//    Output, int &LENGTH, the length of the longest increasing subsequence.
//
//    Output, int SUB[N], contains in entries 1 through LENGTH
//    a longest increasing subsequence of A.
//
{
  int a1;
  int a2;
  int i;
  int j;
  int k;
  int *top;
  int *top_prev;

  if ( n <= 0 )
  {
    length = 0;
    return;
  }

  top = new int[n];
  for ( i = 0; i < n; i++ )
  {
    top[i] = 0;
  }

  top_prev = new int[n];
  for ( i = 0; i < n; i++ )
  {
    top_prev[i] = 0;
  }
  for ( i = 0; i < n; i++ )
  {
    sub[i] = 0;
  }

  length = 0;

  for ( i = 1; i <= n; i++ )
  {
    k = 0;

    for ( j = 1; j <= length; j++ )
    {
      if ( a[i-1] <= a[top[j-1]-1] )
      {
        k = j;
        break;
      }
    }

    if ( k == 0 )
    {
      length = length + 1;
      k = length;
    }

    top[k-1] = i;

    if ( 1 < k )
    {
      top_prev[i-1] = top[k-2];
    }
    else
    {
      top_prev[i-1] = 0;
    }
  }

  j = top[length-1];
  sub[length-1] = a[j-1];

  for ( i = length-1; 1 <= i; i-- )
  {
    j = top_prev[j-1];
    sub[i-1] = a[j-1];
  }

  delete [] top;
  delete [] top_prev;

  return;
}
//****************************************************************************80

int perm_break_count ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_BREAK_COUNT counts the number of "breaks" in a permutation.
//
//  Discussion:
//
//    We begin with a permutation of order N.  We prepend an element
//    labeled "0" and append an element labeled "N+1".  There are now
//    N+1 pairs of neighbors.  A "break" is a pair of neighbors whose
//    value differs by more than 1.  
//
//    The identity permutation has a break count of 0.  The maximum
//    break count is N+1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the permutation.
//
//    Input, int P[N], a permutation, in standard index form.
//
//    Output, int PERM_BREAK_COUNT, the number of breaks in the permutation.
//
{
  bool error;
  int i;
  int value;

  value = 0;
//
//  Make sure the permutation is a legal one.
//  (This is not an efficient way to do so!)
//
  error = perm_check ( n, p );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_BREAK_COUNT - Fatal error!\n";
    cerr << "  The input array does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }

  if ( p[0] != 1 )
  {
    value = value + 1;
  }

  for ( i = 1; i <= n-1; i++ )
  {
    if ( abs ( p[i] - p[i-1] ) != 1 )
    {
      value = value + 1;
    }
  }

  if ( p[n-1] != n )
  {
    value = value + 1;
  }

  return value;
}
//****************************************************************************80

void perm_canon_to_cycle ( int n, int p1[], int p2[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CANON_TO_CYCLE converts a permutation from canonical to cycle form.
//
//  Example:
//
//    Input:
//
//      4 5 2 1 6 3
//
//    Output:
//
//      -4 5 -2 -1 6 3,
//      indicating the cycle structure
//      ( 4, 5 ) ( 2 ) ( 1, 6, 3 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 1, Fundamental Algorithms,
//    Addison Wesley, 1968, page 176.
//
//  Parameters:
//
//    Input, int N, the number of objects permuted.
//
//    Input, int P1[N], the permutation, in canonical form.
//
//    Output, int P2[N], the permutation, in cycle form.
//
{
  int i;
  int pmin;

  for ( i = 0; i < n; i++ )
  {
    p2[i] = p1[i];
  }

  pmin = p2[0] + 1;

  for ( i = 0; i < n; i++ )
  {
    if ( p2[i] < pmin )
    {
      pmin = p2[i];
      p2[i] = -p2[i];
    }
  }
  return;
}
//****************************************************************************80

bool perm_check ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from 1
//    to N occurs among the N entries of the permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the permutation, in standard index form.
//
//    Output, bool PERM_CHECK, is true if the array is NOT a permutation.
//
{
  bool error;
  int ifind;
  int iseek;

  for ( iseek = 1; iseek <= n; iseek++ )
  {
    error = true;

    for ( ifind = 1; ifind <= n; ifind++ )
    {
      if ( p[ifind-1] == iseek )
      {
        error = false;
        break;
      }
    }

    if ( error )
    {
      return true;
    }
  }

  return false;
}
//****************************************************************************80

void perm_cycle ( int n, int p[], int &isgn, int &ncycle, int iopt )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CYCLE analyzes a permutation.
//
//  Discussion:
//
//    The routine will count cycles, find the sign of a permutation,
//    and tag a permutation.
//
//  Example:
//
//    Input:
//
//      N = 9
//      IOPT = 1
//      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
//
//    Output:
//
//      NCYCLE = 3
//      ISGN = +1
//      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N].  On input, P describes a
//    permutation, in the sense that entry I is to be moved to P[I].
//    If IOPT = 0, then P will not be changed by this routine.
//    If IOPT = 1, then on output, P will be "tagged".  That is,
//    one element of every cycle in P will be negated.  In this way,
//    a user can traverse a cycle by starting at any entry I1 of P
//    which is negative, moving to I2 = ABS(P[I1]), then to
//    P[I2], and so on, until returning to I1.
//
//    Output, int &ISGN, the "sign" of the permutation, which is
//    +1 if the permutation is even, -1 if odd.  Every permutation
//    may be produced by a certain number of pairwise switches.
//    If the number of switches is even, the permutation itself is
//    called even.
//
//    Output, int &NCYCLE, the number of cycles in the permutation.
//
//    Input, int IOPT, requests tagging.
//    0, the permutation will not be tagged.
//    1, the permutation will be tagged.
//
{
  bool error;
  int i;
  int i1;
  int i2;
  int is;

  error = perm_check ( n, p );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_CYCLE - Fatal error!\n";
    cerr << "  The input array does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }

  is = 1;
  ncycle = n;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p[i-1];

    while ( i < i1 )
    {
      ncycle = ncycle - 1;
      i2 = p[i1-1];
      p[i1-1] = -i2;
      i1 = i2;
    }

    if ( iopt != 0 )
    {
      is = - i4_sign ( p[i-1] );
    }
    p[i-1] = abs ( p[i-1] ) * i4_sign ( is );
  }

  isgn = 1 - 2 * ( ( n - ncycle ) % 2 );

  return;
}
//****************************************************************************80

void perm_cycle_to_canon ( int n, int p1[], int p2[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CYCLE_TO_CANON converts a permutation from cycle to canonical form.
//
//  Example:
//
//    Input:
//
//      -6 3 1 -5, 4 -2,
//      indicating the cycle structure
//      ( 6, 3, 1 ) ( 5, 4 ) ( 2 )
//
//    Output:
//
//      4 5 2 1 6 3
//
//  Discussion:
//
//    The procedure is to "rotate" the elements of each cycle so that
//    the smallest element is first:
//
//      ( 1, 6, 3 ) ( 4, 5 ) ( 2 )
//
//    and then to sort the cycles in decreasing order of their first
//    (and lowest) element:
//
//      ( 4, 5 ) ( 2 ) ( 1, 6, 3 )
//
//    and then to drop the parentheses:
//
//      4 5 2 1 6 3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 1, Fundamental Algorithms,
//    Addison Wesley, 1968, pages 176.
//
//  Parameters:
//
//    Input, int N, the number of objects permuted.
//
//    Input, int P1[N], the permutation, in cycle form.
//
//    Output, int P2[N], the permutation, in canonical form.
//
{
  int *hi;
  int i;
  int *indx;
  int j;
  int k;
  int *lo;
  int ncycle;
  int next;
  int nhi;
  int nlo;
  int nmin;
  int *pmin;
  int *ptemp;

  hi = new int[n];
  lo = new int[n];
  pmin = new int[n];
  ptemp = new int[n];

  for ( i = 0; i < n; i++ )
  {
    p2[i] = p1[i];
  }
//
//  Work on the next cycle.
//
  nlo = 1;
  ncycle = 0;

  while ( nlo <= n )
//
//  Identify NHI, the last index in this cycle.
//
  {
    ncycle = ncycle + 1;

    nhi = nlo;

    while ( nhi < n )
    {
      if ( p2[nhi] < 0 )
      {
        break;
      }
      nhi = nhi + 1;
    }
//
//  Identify the smallest value in this cycle.
//
    p2[nlo-1] = -p2[nlo-1];
    pmin[ncycle-1] = p2[nlo-1];
    nmin = nlo;

    for ( i = nlo+1; i <= nhi; i++ )
    {
      if ( p2[i-1] < pmin[ncycle-1] )
      {
        pmin[ncycle-1] = p2[i-1];
        nmin = i;
      }
    }
//
//  Rotate the cycle so A_MIN occurs first.
//
    for ( i = nlo; i <= nmin-1; i++ )
    {
      ptemp[i+nhi-nmin] = p2[i-1];
    }
    for ( i = nmin; i <= nhi; i++ )
    {
      ptemp[i+nlo-nmin-1] = p2[i-1];
    }

    lo[ncycle-1] = nlo;
    hi[ncycle-1] = nhi;
//
//  Prepare to operate on the next cycle.
//
    nlo = nhi + 1;
  }
//
//  Compute a sorting index for the cycle minima.
//
  indx = i4vec_sort_heap_index_d ( ncycle, pmin );
//
//  Copy the cycles out of the temporary array in sorted order.
//
  j = 0;

  for ( i = 0; i < ncycle; i++ )
  {
    next = indx[i];
    nlo = lo[next];
    nhi = hi[next];
 
    for ( k = nlo; k <= nhi; k++ )
    {
      j = j + 1;
      p2[j-1] = ptemp[k-1];
    }
  }

  delete [] hi;
  delete [] indx;
  delete [] lo;
  delete [] pmin;
  delete [] ptemp;

  return;
}
//****************************************************************************80

void perm_cycle_to_index ( int n, int p1[], int p2[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CYCLE_TO_INDEX converts a permutation from cycle to standard index form.
//
//  Example:
//
//    Input:
//
//      N = 9
//      P1 = -1, 2, 3, 9, -4, 6, 8, -5, 7
//
//    Output:
//
//      P2 = 2, 3, 9, 6, 7, 8, 5, 4, 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input, int P1[N], the permutation, in cycle form.
//
//    Output, int P2[N], the permutation, in standard index form.
//
{
  int j;
  int k1;
  int k2;
  int k3;

  for ( j = 1; j <= n; j++ )
  {
    k1 = p1[j-1];

    if ( k1 < 0 )
    {
      k1 = -k1;
      k3 = k1;
    }

    if ( j + 1 <= n )
    {
      k2 = p1[j];
      if ( k2 < 0 )
      {
        k2 = k3;
      }
    }
    else
    {
      k2 = k3;
    }

    p2[k1-1] = k2;

  }

  return;
}
//****************************************************************************80

int perm_distance ( int n, int a[], int b[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_DISTANCE computes the Ulam metric distance of two permutations.
//
//  Discussion:
//
//    If we let N be the order of the permutations A and B, and L(P) be
//    the length of the longest ascending subsequence of a permutation P,
//    then the Ulam metric distance between A and B is
//
//      N - L ( A * inverse ( B ) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the permutation.
//
//    Input, int A[N], B[N], the permutations to be examined.
//
//    Output, int K, the Ulam metric distance between A and B.
//
{
  int *binv;
  int *c;
  int i;
  int length;
  int *sub;

  binv = new int[n];
  c = new int[n];
  sub = new int[n];

  for ( i = 0; i < n; i++ )
  {
    binv[i] = b[i];
  }

  perm_inverse ( n, binv );

  perm_mul ( n, a, binv, c );

  perm_ascend ( n, c, length, sub );

  delete [] binv;
  delete [] c;
  delete [] sub;

  return ( n - length );
}
//****************************************************************************80

int perm_fixed_enum ( int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_FIXED_ENUM enumerates the permutations of N objects with M fixed.
//
//  Discussion:
//
//    A permutation of N objects with M fixed is a permutation in which
//    exactly M of the objects retain their original positions.  If
//    M = 0, the permutation is a "derangement".  If M = N, the
//    permutation is the identity.
//
//    The formula is:
//
//      F(N,M) = ( N! / M! ) * ( 1 - 1/1! + 1/2! - 1/3! ... 1/(N-M)! )
//             = COMB(N,M) * D(N-M)
//
//    where D(N-M) is the number of derangements of N-M objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//    N should be at least 1.
//
//    Input, int M, the number of objects that retain their
//    position.  M should be between 0 and N.
//
//    Output, int PERM_FIXED_ENUM, the number of derangements of N objects
//    in which M objects retain their positions.
//
{
  int fnm;

  if ( n <= 0 )
  {
    fnm = 1;
  }
  else if ( m < 0 )
  {
    fnm = 0;
  }
  else if ( n < m )
  {
    fnm = 0;
  }
  else if ( m == n )
  {
    fnm = 1;
  }
  else if ( n == 1 )
  {
    if ( m == 1 )
    {
      fnm = 1;
    }
    else
    {
      fnm = 0;
    }
  }
  else
  {
    fnm = i4_choose ( n, m ) * derange_enum ( n - m );
  }

  return fnm;
}
//****************************************************************************80

void perm_free ( int npart, int ipart[], int nfree, int ifree[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_FREE reports the unused items in a partial permutation.
//
//  Discussion:
//
//    It is assumed that the N objects being permuted are the integers
//    from 1 to N, and that IPART contains a "partial" permutation, that
//    is, the NPART entries of IPART represent the beginning of a
//    permutation of all N items.
//
//    The routine returns in IFREE the items that have not been used yet.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NPART, the number of entries in IPART.  NPART may be 0.
//
//    Input, int IPART[NPART], the partial permutation, which should
//    contain, at most once, some of the integers between 1 and
//    NPART+NFREE.
//
//    Input, int NFREE, the number of integers that have not been
//    used in IPART.  This is simply N - NPART.  NFREE may be zero.
//
//    Output, int IFREE[NFREE], the integers between 1 and NPART+NFREE
//    that were not used in IPART.
//
{
  int i;
  int j;
  int k;
  int match;
  int n;

  n = npart + nfree;

  if ( npart < 0 )
  {
    cerr << "\n";
    cerr << "PERM_FREE - Fatal error!\n";
    cerr << "  NPART < 0.\n";
    exit ( 1 );
  }
  else if ( npart == 0 )
  {
    i4vec_indicator ( n, ifree );
  }
  else if ( nfree < 0 )
  {
    cerr << "\n";
    cerr << "PERM_FREE - Fatal error!\n";
    cerr << "  NFREE < 0.\n";
    exit ( 1 );
  }
  else if ( nfree == 0 )
  {
    return;
  }
  else
  {
    k = 0;

    for ( i = 1; i <= n; i++ )
    {
      match = -1;

      for ( j = 0; j < npart; j++ )
      {
        if ( ipart[j] == i )
        {
          match = j;
          break;
        }
      }

      if ( match == -1 )
      {
        k = k + 1;

        if ( nfree < k )
        {
          cerr << "\n";
          cerr << "PERM_FREE - Fatal error!\n";
          cerr << "  The partial permutation is illegal.\n";
          cerr << "  Technically, because NFREE < K.\n";
          cerr << "  N     = " << n << "\n";
          cerr << "  NPART = " << npart << "\n";
          cerr << "  NFREE = " << nfree << "\n";
          cerr << "  K =     " << k << "\n";
          cerr << "\n";
          cerr << "  The partial permutation:\n";
          cerr << "\n";
          for ( i = 0; i < npart; i++ )
          {
            cerr << setw(2) << ipart[i] << "  ";
          }
          cerr << "\n";
          exit ( 1 );
        }
        ifree[k-1] = i;
      }
    }
  }

  return;
}
//****************************************************************************80

void perm_index_to_cycle ( int n, int p1[], int p2[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INDEX_TO_CYCLE converts a permutation from standard index to cycle form.
//
//  Example:
//
//    Input:
//
//      N = 9
//      P1 = 2, 3, 9, 6, 7, 8, 5, 4, 1
//
//    Output:
//
//      P2 = -1, 2, 3, 9, -4, 6, 8, -5, 7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input, int P1[N], the permutation, in standard index form.
//
//    Output, int P2[N], the permutation, in cycle form.
//
{
  int i;
  int j;
  int k;

  i = 0;
  j = 1;

  while ( j <= n )
  {
    if ( p1[j-1] < 0 )
    {
      j = j + 1;
    }
    else
    {
      k = j;

      i = i + 1;
      p2[i-1] = - k;

      while ( p1[k-1] != j )
      {
        i = i + 1;
        p2[i-1] = p1[k-1];
        p1[k-1] = - p1[k-1];
        k = abs ( p1[k-1] );
      }

      p1[k-1] = - p1[k-1];
    }
  }

  for ( i = 0; i < n; i++ )
  {
    p1[i] = abs ( p1[i] );
  }

  return;
}
//****************************************************************************80

void perm_ins ( int n, int p[], int ins[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INS computes the inversion sequence of a permutation.
//
//  Discussion:
//
//    For a given permutation P acting on objects 1 through N, the inversion
//    sequence INS is defined as:
//
//      INS(1) = 0
//      INS(I) = number of values J < I for which P(I) < P(J).
//
//  Example:
//
//    Input:
//
//      ( 3, 5, 1, 4, 2 )
//
//    Output:
//
//      ( 0, 0, 2, 1, 3 )
//
//    The original permutation can be recovered from the inversion sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input, int P[N], the permutation, in standard index form.
//    The I-th item has been mapped to P(I).
//
//    Output, int INS[N], the inversion sequence of the permutation.
//
{
  bool error;
  int i;
  int j;

  error = perm_check ( n, p );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_INS - Fatal error!\n";
    cerr << "  The input array does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    ins[i] = 0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      if ( p[i] < p[j] )
      {
        ins[i] = ins[i] + 1;
      }
    }
  }

  return;
}
//****************************************************************************80

void perm_inverse ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE inverts a permutation "in place".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    On output, P describes the inverse permutation
//
{
  bool error;
  int i;
  int i0;
  int i1;
  int i2;
  int is;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "PERM_INVERSE - Fatal error!\n";
    cerr << "  Input value of N = " << n << "\n";
    exit ( 1 );
  }

  error = perm_check ( n, p );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_INVERSE - Fatal error!\n";
    cerr << "  The input array does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }

  is = 1;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p[i-1];

    while ( i < i1 )
    {
      i2 = p[i1-1];
      p[i1-1] = -i2;
      i1 = i2;
    }

    is = -i4_sign ( p[i-1] );
    p[i-1] = abs ( p[i-1] ) * i4_sign ( is );

  }

  for ( i = 1; i <= n; i++ )
  {
    i1 = -p[i-1];

    if ( 0 <= i1 )
    {
      i0 = i;

      for ( ; ; )
      {
        i2 = p[i1-1];
        p[i1-1] = i0;

        if ( i2 < 0 )
        {
          break;
        }

        i0 = i1;
        i1 = i2;
      }
    }
  }

  return;
}
//****************************************************************************80

void perm_inverse2 ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE2 inverts a permutation "in place".
//
//  Discussion:
//
//    The routine needs no extra vector storage in order to compute the
//    inverse of a permutation.
//
//    This feature might be useful if the permutation is large.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of objects in the permutation.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    On output, the inverse permutation.
//
{
  bool error;
  int i;
  int ii;
  int j;
  int k;
  int m;

  error = perm_check ( n, p );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_INVERSE2 - Fatal error!\n";
    cerr << "  The input array does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }

  for ( ii = 1; ii <= n; ii++ )
  {
    m = n + 1 - ii;
    i = p[m-1];

    if ( i < 0 )
    {
      p[m-1] = -i;
    }
    else if ( i != m )
    {
      k = m;

      for ( ; ; )
      {
        j = p[i-1];
        p[i-1] = -k;

        if ( j == m )
        {
          p[m-1] = i;
          break;
        }

        k = i;
        i = j;
      }
    }
  }

  return;
}
//****************************************************************************80

int *perm_inverse3_new ( int n, int perm[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE3 produces the inverse of a given permutation.
//
//  Discussion:
//
//    The input permutation is assumed to be 0-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items permuted.
//
//    Input, int PERM[N], a permutation.
//
//    Output, int PERM_INVERSE3_NEW[N], the inverse permutation.
//
{
  int i;
  int *perm_inv;

  perm_inv = new int[n];

  for ( i = 0; i < n; i++ )
  {
    perm_inv[perm[i]] = i;
  }

  return perm_inv;
}
//****************************************************************************80

void perm_lex_next ( int n, int p[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_LEX_NEXT generates permutations in lexical order, one at a time.
//
//  Example:
//
//    N = 3
//
//    1   1 2 3
//    2   1 3 2
//    3   2 1 3
//    4   2 3 1
//    5   3 1 2
//    6   3 2 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Mok-Kong Shen.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Mok-Kong Shen,
//    Algorithm 202: Generation of Permutations in Lexicographical Order,
//    Communications of the ACM,
//    Volume 6, September 1963, page 517.
//
//  Parameters:
//
//    Input, int N, the number of elements being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    The user should not alter the elements of Pbetween successive
//    calls.  The routine uses the input value of P to determine
//    the output value.
//
//    Input/output, bool &MORE.
//    On the first call, the user should set MORE =FALSE which signals
//    the routine to do initialization.
//    On return, if MORE is TRUE then another permutation has been
//    computed and returned, while if MORE is FALSE there are no more
//    permutations.
//
{
  int i;
  int j;
  int k;
  int temp;
  int u;
  int w;
//
//  Initialization.
//
  if ( !more )
  {
    i4vec_indicator ( n, p );
    more = true;
  }
  else
  {
    if ( n <= 1 )
    {
      more = false;
      return;
    }

    w = n;

    while ( p[w-1] < p[w-2] )
    {
      if ( w == 2 )
      {
        more = false;
        return;
      }

      w = w - 1;
    }

    u = p[w-2];

    for ( j = n; w <= j; j-- )
    {
      if ( u < p[j-1] )
      {
        p[w-2] = p[j-1];
        p[j-1] = u;

        for ( k = 0; k <= ( n - w - 1 ) / 2; k++ )
        {
          temp = p[n-k-1];
          p[n-k-1] = p[w+k-1];
          p[w+k-1] = temp;
        }
        return;
      }
    }
  }

  return;
}
//****************************************************************************80

void perm_mul ( int n, int p1[], int p2[], int p3[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_MUL "multiplies" two permutations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the permutations.
//
//    Input, int P1[N], P2[N], the permutations, in standard index form.
//
//    Output, int P3[N], the product permutation.
//
{
  bool error;
  int i;

  error = perm_check ( n, p1 );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_MUL - Fatal error!\n";
    cerr << "  The input array P1 does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }

  error = perm_check ( n, p2 );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_MUL - Fatal error!\n";
    cerr << "  The input array P2 does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    p3[i] = p2[p1[i]-1];
  }

  return;
}
//****************************************************************************80

void perm_next ( int n, int p[], bool &more, bool &even )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_NEXT computes all of the permutations of N objects, one at a time.
//
//  Discussion:
//
//    The routine is initialized by calling with MORE = TRUE, in which case
//    it returns the identity permutation.
//
//    If the routine is called with MORE = FALSE, then the successor of the
//    input permutation is computed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    On the first call, the input value is unimportant.
//    On subsequent calls, the input value should be the same
//    as the output value from the previous call.  In other words, the
//    user should just leave P alone.
//    On output, contains the "next" permutation.
//
//    Input/output, bool &MORE.
//    Set MORE = FALSE before the first call.
//    MORE will be reset to TRUE and a permutation will be returned.
//    Each new call produces a new permutation until
//    MORE is returned FALSE.
//
//    Input/output, bool &EVEN.
//    The input value of EVEN should simply be its output value from the
//    previous call; (the input value on the first call doesn't matter.)
//    On output, EVEN is TRUE if the output permutation is even, that is,
//    involves an even number of transpositions.
//
{
  int i;
  int i1;
  int ia;
  int id;
  int is;
  int j;
  int l;
  int m;

  if ( !more )
  {
    i4vec_indicator ( n, p );

    more = true;
    even = true;

    if ( n == 1 )
    {
      more = false;
      return;
    }

    if ( p[n-1] != 1 || p[0] != 2 + ( n % 2 ) )
    {
      return;
    }

    for ( i = 1; i <= n-3; i++ )
    {
      if ( p[i] != p[i-1] + 1 )
      {
        return;
      }
    }

    more = false;
  }
  else
  {
    if ( n == 1 )
    {
      p[0] = 0;
      more = false;
      return;
    }

    if ( even )
    {
      ia = p[0];
      p[0] = p[1];
      p[1] = ia;
      even = false;

      if ( p[n-1] != 1 || p[0] != 2 + ( n % 2 ) )
      {
        return;
      }

      for ( i = 1; i <= n-3; i++ )
      {
        if ( p[i] != p[i-1] + 1 )
        {
          return;
        }
      }

      more = false;
      return;
    }
    else
    {
      more = false;

      is = 0;

      for ( i1 = 2; i1 <= n; i1++ )
      {
        ia = p[i1-1];
        i = i1 - 1;
        id = 0;

        for ( j = 1; j <= i; j++ )
        {
          if ( ia < p[j-1] )
          {
            id = id + 1;
          }
        }

        is = id + is;

        if ( id != i * ( is % 2 ) )
        {
          more = true;
          break;
        }
      }

      if ( !more )
      {
        p[0] = 0;
        return;
      }
    }

    m = ( ( is + 1 ) % 2 ) * ( n + 1 );

    for ( j = 1; j <= i; j++ )
    {
      if ( i4_sign ( p[j-1] - ia ) != i4_sign ( p[j-1] - m ) )
      {
        m = p[j-1];
        l = j;
      }
    }

    p[l-1] = ia;
    p[i1-1] = m;
    even = true;
  }

  return;
}
//****************************************************************************80

void perm_next2 ( int n, int p[], bool &done )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_NEXT2 generates all the permutations of N objects.
//
//  Discussion:
//
//    The routine generates the permutations one at a time.  It uses a
//    particular ordering of permutations, generating them from the first
//    (which is the identity permutation) to the N!-th.  The same ordering
//    is used by the routines PERM_RANK and PERM_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 October 2006
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input, int N, the number of elements in the set to be permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//
//    Input/output, bool &DONE.  The user should set the input value of
//    DONE only once, before the first call to compute the permutations.
//    The user should set DONE to TRUE, which signals the routine
//    that it is to initialize itself.
//    Thereafter, the routine will set DONE to FALSE and will
//    compute a new permutation on each call.
//    However, when there are no more permutations to compute, the
//    routine will not return a new permutation, but instead will
//    return DONE with the value TRUE.  At this point, all the
//    permutations have been computed.
//
{
  static int *active = NULL;
  int i;
  static int *idir = NULL;
  static int *invers = NULL;
  int j;
  int nactiv;
//
//  An input value of FALSE for DONE is assumed to mean a new
//  computation is beginning.
//
  if ( done )
  {
    i4vec_indicator ( n, p );

    if ( active )
    {
      delete [] active;
    }
    active = new int[n];
    if ( idir )
    {
      delete [] idir;
    }
    idir = new int[n];
    if ( invers )
    {
      delete [] invers;
    }
    invers = new int[n];

    for ( i = 0; i < n; i++ )
    {
      invers[i] = p[i];
    }
    for ( i = 0; i < n; i++ )
    {
      idir[i] = -1;
    }
    active[0] = 0;
    for ( i = 1; i < n; i++ )
    {
      active[i] = 1;
    }
//
//  Set the DONE flag to FALSE, signifying there are more permutations
//  to come.  Except, of course, that we must take care of the trivial case!
//
    if ( 1 < n )
    {
      done = false;
    }
    else
    {
      done = true;
    }
  }
//
//  Otherwise, assume we are in a continuing computation
//
  else
  {
    nactiv = 0;

    for ( i = 1; i <= n; i++ )
    {
      if ( active[i-1] != 0 )
      {
        nactiv = i;
      }
    }

    if ( nactiv <= 0 )
    {
      done = true;
    }
    else
    {
      j = invers[nactiv-1];

      p[j-1] = p[j+idir[nactiv-1]-1];
      p[j+idir[nactiv-1]-1] = nactiv;

      invers[nactiv-1] = invers[nactiv-1] + idir[nactiv-1];
      invers[p[j-1]-1] = j;

      if ( j + 2 * idir[nactiv-1] < 1 || n < j + 2 * idir[nactiv-1] )
      {
        idir[nactiv-1] = - idir[nactiv-1];
        active[nactiv-1] = 0;
      }
      else if ( nactiv < p[j+2*idir[nactiv-1]-1] )
      {
        idir[nactiv-1] = - idir[nactiv-1];
        active[nactiv-1] = 0;
      }

      for ( i = nactiv; i < n; i++ )
      {
        active[i] = 1;
      }
    }
  }

  if ( done )
  {
    delete [] active;
    active = NULL;
    delete [] idir;
    idir = NULL;
    delete [] invers;
    invers = NULL;
  }

  return;
}
//****************************************************************************80

void perm_next3 ( int n, int p[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_NEXT3 computes all of the permutations of N objects, one at a time.
//
//  Discussion:
//
//    The routine is initialized by calling with MORE = TRUE in which case
//    it returns the identity permutation.
//
//    If the routine is called with MORE = FALSE then the successor of the
//    input permutation is computed.
//
//    Trotter's algorithm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 July 2003
//
//  Author:
//
//    Original FORTRAN77 version by H Trotter,
//    C++ version by John Burkardt.
//
//  Reference:
//
//    H Trotter,
//    PERM, Algorithm 115,
//    Communications of the Association for Computing Machinery,
//    Volume 5, 1962, pages 434-435.
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    If MORE is TRUE, then P is assumed to contain the
//    "previous" permutation, and on P(I) is the value
//    of the I-th object under the next permutation.
//    Otherwise, P will be set to the "first" permutation.
//
//    Input/output, bool &MORE.
//    Set MORE = FALSE before first calling this routine.
//    MORE will be reset to TRUE and a permutation will be returned.
//    Each new call produces a new permutation until MORE is returned FALSE
//
{
  int i;
  int m2;
  int n2;
  static int nfact = 0;
  int q;
  static int rank = 0;
  int s;
  int t;
  int temp;

  if ( !more )
  {
    i4vec_indicator ( n, p );
    more = true;
    rank = 1;

    nfact = i4_factorial ( n );
  }
  else
  {
    n2 = n;
    m2 = rank;
    s = n;

    for ( ; ; )
    {
      q = m2 % n2;
      t = m2 % ( 2 * n2 );

      if ( q != 0 )
      {
        break;
      }

      if ( t == 0 )
      {
        s = s - 1;
      }

      m2 = m2 / n2;
      n2 = n2 - 1;

    }

    if ( q == t )
    {
      s = s - q;
    }
    else
    {
      s = s + q - n2;
    }

    temp = p[s-1];
    p[s-1] = p[s];
    p[s] = temp;

    rank = rank + 1;

    if ( rank == nfact )
    {
      more = false;
    }
  }

  return;
}
//****************************************************************************80

void perm_print ( int n, int p[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_PRINT prints a permutation.
//
//  Discussion:
//
//    The permutation is assumed to be zero-based.
//
//  Example:
//
//    Input:
//
//      P = 6 1 2 0 4 2 5
//
//    Printed output:
//
//      "This is the permutation:"
//
//      0 1 2 3 4 5 6
//      6 1 2 0 4 2 5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects permuted.
//
//    Input, int P[N], the permutation, in standard index form.
//
//    Input, string TITLE, a title.
//    If no title is supplied, then only the permutation is printed.
//
{
  int i;
  int ihi;
  int ilo;
  int inc = 20;

  if ( s_len_trim ( title ) != 0 )
  {
    cout << "\n";
    cout << title << "\n";

    for ( ilo = 0; ilo < n; ilo = ilo + inc )
    {
      ihi = ilo + inc;
      if ( n < ihi ) 
      {
        ihi = n;
      }
      cout << "\n";
      cout << "  ";
      for ( i = ilo; i < ihi; i++ )
      {
        cout << setw(4) << i;
      }
      cout << "\n";
      cout << "  ";
      for ( i = ilo; i < ihi; i++ )
      {
        cout << setw(4) << p[i];
      }
      cout << "\n";
    }
  }
  else
  {
    for ( ilo = 0; ilo < n; ilo = ilo + inc )
    {
      ihi = ilo + inc;
      if ( n < ihi ) 
      {
        ihi = n;
      }
      cout << "  ";
      for ( i = ilo; i < ihi; i++ )
      {
        cout << setw(4) << p[i];
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

void perm_random ( int n, int &seed, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_RANDOM selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int P[N], a permutation of (1, 2, ..., N).
//
{
  int i;
  int j;
  int k;
 
  for ( i = 0; i < n; i++ )
  {
    p[i] = i + 1;
  }

  for ( i = 1; i <= n; i++ )
  {
    j = i4_uniform ( i, n, seed );
    k = p[i-1];
    p[i-1] = p[j-1];
    p[j-1] = k;
  }
 
  return;
}
//****************************************************************************80

void perm_random2 ( int n, int &seed, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_RANDOM2 selects a random permutation of N objects.
//
//  Discussion:
//
//    The input values of P are used as labels; that is, the I-th
//    object is labeled P[I-1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Input/utput, int P[N], a permutation of the original values in P.
//
{
  int i;
  int j;
  int k;

  for ( i = 1; i <= n; i++ )
  {
    j = i4_uniform ( i, n, seed );
    k = p[i-1];
    p[i-1] = p[j-1];
    p[j-1] = k;
  }
 
  return;
}
//****************************************************************************80

void perm_random3 ( int n, int &seed, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_RANDOM3 selects a random permutation of N elements.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by James Filliben.
//    C++ version by John Burkardt
//
//  Reference:
//
//    K L Hoffman, D R Shier,
//    Algorithm 564,
//    A Test Problem Generator for Discrete Linear L1 Approximation Problems,
//    ACM Transactions on Mathematical Software,
//    Volume 6, Number 4, December 1980, pages 615-617.
//
//  Parameters:
//
//    Input, int N, the number of elements of the array.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int P[N], a permutation, in standard index form.
//
{
  int i;
  int j;
  int temp;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "PERM_RANDOM3 - Fatal error!\n";
    cerr << "  Illegal input value of N  = " << n << "\n";
    cerr << "  N must be at least 1!\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    p[0] = 1;
    return;
  }

  i4vec_indicator ( n, p );

  for ( i = 1; i <= n; i++ )
  {
    j = i + i4_uniform ( 1, n, seed );

    if ( n < j )
    {
      j = j - n;
    }

    if ( i != j )
    {
      temp = p[j-1];
      p[j-1] = p[i-1];
      p[i-1] = temp;
    }
  }

  return;
}
//****************************************************************************80

int perm_rank ( int n, int p[], int invers[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_RANK computes the rank of a given permutation.
//
//  Discussion:
//
//    This is the same as asking for the step at which PERM_NEXT2
//    would compute the permutation.  The value of the rank will be
//    between 1 and N!.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input, int N, the number of elements in the set that
//    is permuted by P.
//
//    Input, int P[N], a permutation, in standard index form.
//
//    Output, int INVERS[N], the inverse permutation of P.
//    It is computed as part of the algorithm, and may be of use
//    to the user.  INVERS(P(I)) = I for each entry I.
//
//    Output, int PERM_RANK, the rank of the permutation.  This
//    gives the order of the given permutation in the set of all
//    the permutations on N elements.
//
{
  int count;
  bool error;
  int i;
  int j;
  int rank;
  int rem;
//
//  Make sure the permutation is a legal one.
//  (This is not an efficient way to do so!)
//
  error = perm_check ( n, p );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_RANK - Fatal error!\n";
    cerr << "  The input array does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }
//
//  Compute the inverse permutation.
//
  for ( i = 0; i < n; i++ )
  {
    invers[i] = p[i];
  }

  perm_inverse2 ( n, invers );

  rank = 0;

  for ( i = 1; i <= n; i++)
  {

    count = 0;

    for ( j = 0; j < invers[i-1]; j++ )
    {
      if ( p[j] < i )
      {
        count = count + 1;
      }
    }

    if ( ( rank % 2 ) == 1 )
    {
      rem = count;
    }
    else
    {
      rem = i - 1 - count;
    }

    rank = i * rank + rem;
  }

  rank = rank + 1;

  return rank;
}
//****************************************************************************80

int perm_sign ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_SIGN returns the sign of a permutation.
//
//  Discussion:
//
//    A permutation can always be replaced by a sequence of pairwise
//    transpositions.  A given permutation can be represented by
//    many different such transposition sequences, but the number of
//    such transpositions will always be odd or always be even.
//    If the number of transpositions is even or odd, the permutation is
//    said to be even or odd.
//
//  Example:
//
//    Input:
//
//      N = 9
//      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
//
//    Output:
//
//      PERM_SIGN = +1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2012
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of objects permuted.
//
//    Input, int P[N], a permutation, in standard index form.
//
//    Output, int PERM_SIGN, the "sign" of the permutation.
//    +1, the permutation is even,
//    -1, the permutation is odd.
//
{
  bool error;
  int i;
  int j;
  int p_sign;
  int *q;
  int temp;

  error = perm_check ( n, p );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_SIGN - Fatal error!\n";
    cerr << "  The input array does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }
//
//  Make a temporary copy of the permutation.
//
  q = new int[n];
  for ( i = 0; i < n; i++ )
  {
    q[i] = p[i];
  }
//
//  Start with P_SIGN indicating an even permutation.
//  Restore each element of the permutation to its correct position,
//  updating P_SIGN as you go.
//
  p_sign = 1;

  for ( i = 1; i <= n - 1; i++ )
  {
    j = i4vec_index ( n, q, i );

    if ( j != i - 1 )
    {
      temp = q[i-1];
      q[i-1] = q[j];
      q[j] = temp;

      p_sign = - p_sign;
    }
  }
  
  delete [] q;
  
  return p_sign;
}
//****************************************************************************80

void perm_to_equiv ( int n, int p[], int &npart, int jarray[], int iarray[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_TO_EQUIV computes the partition induced by a permutation.
//
//  Example:
//
//    Input:
//
//      N = 9
//      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
//
//    Output:
//
//      NPART = 3
//      JARRAY = 4, 3, 2
//      IARRAY = 1, 1, 1, 2  3  2  3  2, 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input, int P[N], a permutation, in standard index form.
//
//    Output, int &NPART, number of subsets in the partition.
//
//    Output, int JARRAY[N].  JARRAY(I) is the number of elements
//    in the I-th subset of the partition.
//
//    Output, int IARRAY[N].  IARRAY(I) is the class to which
//    element I belongs.
//
{
  bool error;
  int i;
  int j;
  int k;

  error = perm_check ( n, p );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_TO_EQUIV - Fatal error!\n";
    cerr << "  The input array does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }
//
//  Initialize.
//
  for ( i = 0; i < n; i++ )
  {
    iarray[i] = 0;
  }
  
  for ( i = 0; i < n; i++ )
  {
    jarray[i] = 0;
  }

  npart = 0;
//
//  Search for the next item J which has not been assigned a subset/orbit.
//
  for ( j = 1; j <= n; j++ )
  {
    if ( iarray[j-1] != 0 )
    {
      continue;
    }
//
//  Begin a new subset/orbit.
//
    npart = npart + 1;
    k = j;
//
//  Add the item to the subset/orbit.
//
    for ( ; ; )
    {
      jarray[npart-1] = jarray[npart-1] + 1;
      iarray[k-1] = npart;
//
//  Apply the permutation.  If the permuted object isn't already in the
//  subset/orbit, add it.
//
      k = p[k-1];

      if ( iarray[k-1] != 0 )
      {
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

void perm_to_ytb ( int n, int p[], int lambda[], int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_TO_YTB converts a permutation to a Young tableau.
//
//  Discussion:
//
//    The mapping is not invertible.  In most cases, several permutations
//    correspond to the same tableau.
//
//  Example:
//
//    N = 7
//    P = 7 2 4 1 5 3 6
//
//    YTB =
//      1 2 3 6
//      4 5
//      7
//
//    LAMBDA = 4 2 1 0 0 0 0
//
//    A = 1 1 1 2 2 1 3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer to be partitioned.
//
//    Input, int P[N], a permutation, in standard index form.
//
//    Output, int LAMBDA[N].  LAMBDA[I] is the length of the I-th row.
//
//    Output, int A[N].  A[I] is the row containing I.
//
{
  bool another;
  int compare;
  int i;
  int put_index;
  int put_row;
  int put_value;
//
//  Initialize.
//
  for ( i = 0; i < n; i++ )
  {
    lambda[i] = 0;
  }
  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
//
//  Now insert each item of the permutation.
//
  for ( put_index = 1; put_index <= n; put_index++ )
  {
    put_value = p[put_index-1];
    put_row = 1;

    for ( ; ; )
    {
      another = false;

      for ( compare = put_value+1; compare <= n; compare++ )
      {
        if ( a[compare-1] == put_row )
        {
          another = true;
          a[put_value-1] = put_row;
          a[compare-1] = 0;
          put_value = compare;
          put_row = put_row + 1;
          break;
        }
      }

      if ( !another )
      {
        break;
      }
    }

    a[put_value-1] = put_row;
    lambda[put_row-1] = lambda[put_row-1] + 1;

  }

  return;
}
//****************************************************************************80

void perm_unrank ( int n, int rank, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNRANK "unranks" a permutation.
//
//  Discussion:
//
//    That is, given a rank, it computes the corresponding permutation.
//    This is the same as asking for the permutation which PERM_NEXT2
//    would compute at the RANK-th step.
//
//    The value of the rank should be between 1 and N!.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2004
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input, int N, the number of elements in the set.
//
//    Input, int RANK, the desired rank of the permutation.  This
//    gives the order of the given permutation in the set of all
//    the permutations on N elements, using the ordering of PERM_NEXT2.
//
//    Output, int P[N], the permutation, in standard index form.
//
{
  int i;
  int icount;
  int iprev;
  int irem;
  int j;
  int jdir;
  int jrank;

  for ( i = 0; i < n; i++ )
  {
    p[i] = 0;
  }

  if ( rank < 1 || i4_factorial ( n ) < rank )
  {
    cerr << "\n";
    cerr << "PERM_UNRANK - Fatal error!\n";
    cerr << "  Illegal input value for RANK.\n";
    cerr << "  RANK must be between 1 and N!,\n";
    cerr << "  but the input value is " << rank << "\n";
    exit ( 1 );
  }

  jrank = rank - 1;

  for ( i = 1; i <= n; i++ )
  {
    iprev = n + 1 - i;
    irem = jrank % iprev;
    jrank = jrank / iprev;

    if ( ( jrank % 2 ) == 1 )
    {
      j = 0;
      jdir = 1;
    }
    else
    {
      j = n + 1;
      jdir = -1;
    }

    icount = 0;

    for ( ; ; )
    {
      j = j + jdir;

      if ( p[j-1] == 0 )
      {
        icount = icount + 1;
      }

      if ( irem < icount )
      {
        break;
      }
    }
    p[j-1] = iprev;
  }

  return;
}
//****************************************************************************80

void perrin ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERRIN returns the first N values of the Perrin sequence.
//
//  Discussion:
//
//    The Perrin sequence has the initial values:
//
//      P(0) = 3
//      P(1) = 0
//      P(2) = 2
//
//    and subsequent entries are generated by the recurrence
//
//      P(I+1) = P(I-1) + P(I-2)
//
//    Note that if N is a prime, then N must evenly divide P(N).
//
//  Example:
//
//    0   3
//    1   0
//    2   2
//    3   3
//    4   2
//    5   5
//    6   5
//    7   7
//    8  10
//    9  12
//   10  17
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Stewart,
//    "A Neglected Number",
//    Scientific American, Volume 274, pages 102-102, June 1996.
//
//    Ian Stewart,
//    Math Hysteria,
//    Oxford, 2004.
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 1999.
//
//  Parameters:
//
//    Input, integer N, the number of terms.
//
//    Output, integer P(N), the terms 0 through N-1 of the sequence.
//
{
  int i;

  if ( n < 1 )
  {
    return;
  }

  p[0] = 3;

  if ( n < 2 )
  {
    return;
  }

  p[1] = 0;

  if ( n < 3 )
  {
    return;
  }
 
  p[2] = 2;

  for ( i = 4; i <= n; i++ )
  {
    p[i-1] = p[i-3] + p[i-4];
  }

  return;
}
//****************************************************************************80

bool pord_check ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    PORD_CHECK checks a matrix representing a partial ordering.
//
//  Discussion:
//
//    The array A is supposed to represent a partial ordering of
//    the elements of a set of N objects.
//
//    For distinct indices I and J, the value of A(I,J) is:
//
//      1, if I << J
//      0, otherwise ( I and J may be unrelated, or perhaps J << I).
//
//    Diagonal elements of A are ignored.
//
//    This routine checks that the values of A do represent
//    a partial ordering.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements in the set.
//
//    Input, int A[N*N], the partial ordering.  A[I+J*N] is
//    1 if I is less than J in the partial ordering,
//    0 otherwise.
//
//    Output, bool PORD_CHECK, is true if an error was detected.
//
{
  int i;
  int j;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "PORD_CHECK - Fatal error!\n";
    cerr << "  N must be positive, but N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      if ( 0 < a[i+j*n] )
      {
        if ( 0 < a[j+i*n] )
        {
          cerr << "\n";
          cerr << "PORD_CHECK - Fatal error!\n";
          cerr << "  For indices I = " << i << "\n";
          cerr << "  and J = " << j << "\n";
          cerr << "  A(I,J) = " << a[i+j*n] << "\n";
          cerr << "  A(J,I) = " << a[j+i*n] << "\n";
          exit ( 1 );
        }
      }
    }
  }

  return false;
}
//****************************************************************************80

int power_mod ( int a, int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_MOD computes mod ( A^N, M ).
//
//  Discussion:
//
//    Some programming tricks are used to speed up the computation, and to
//    allow computations in which A**N is much too large to store in a
//    real word.
//
//    First, for efficiency, the power A**N is computed by determining
//    the binary expansion of N, then computing A, A**2, A**4, and so on
//    by repeated squaring, and multiplying only those factors that
//    contribute to A^N.
//
//    Secondly, the intermediate products are immediately "mod'ed", which
//    keeps them small.
//
//    For instance, to compute mod ( A^13, 11 ), we essentially compute
//
//       13 = 1 + 4 + 8
//
//       A**13 = A * A^4 * A^8
//
//       mod ( A^13, 11 ) = mod ( A, 11 ) * mod ( A^4, 11 ) * mod ( A^8, 11 ).
//
//    Fermat's little theorem says that if P is prime, and A is not divisible
//    by P, then ( A^(P-1) - 1 ) is divisible by P.
//
//  Example:
//
//     A  N  M  X
//
//    13  0 31  1
//    13  1 31 13
//    13  2 31 14
//    13  3 31 27
//    13  4 31 10
//    13  5 31  6
//    13  6 31 16
//    13  7 31 22
//    13  8 31  7 
//    13  9 31 29
//    13 10 31  5
//    13 11 31  3
//    13 12 31  8
//    13 13 31 11
//    13 14 31 19
//    13 15 31 30
//    13 16 31 18
//    13 17 31 17
//    13 18 31  4
//    13 19 31 21
//    13 20 31 25
//    13 21 31 15
//    13 22 31  9
//    13 23 31 24
//    13 24 31  2
//    13 25 31 26
//    13 26 31 28
//    13 27 31 23
//    13 28 31 20
//    13 29 31 12
//    13 30 31  1
//    13 31 31 13
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the base of the expression to be tested.
//    A should be nonnegative.
//
//    Input, int N, the power to which the base is raised.
//    N should be nonnegative.
//
//    Input, int M, the divisor against which the expression is tested.
//    M should be positive.
//
//    Output, int POWER_MOD, the remainder when A^N is divided by M.
//
{
  long long int a_square2;
  int d;
  long long int m2;
  int x;
  long long int x2;

  if ( a < 0 )
  {
    return -1;
  }

  if ( m <= 0 )
  {
    return -1;
  }

  if ( n < 0 )
  {
    return -1;
  }
//
//  A_SQUARE2 contains the successive squares of A.
//
  a_square2 = ( long long int ) a;
  m2 = ( long long int ) m;
  x2 = ( long long int ) 1;

  while ( 0 < n )
  {
    d = n % 2;

    if ( d == 1 )
    {
      x2 = ( x2 * a_square2 ) % m2;
    }

    a_square2 = ( a_square2 * a_square2 ) % m2;
    n = ( n - d ) / 2;
  }
//
//  Ensure that 0 <= X.
//
  while ( x2 < 0 )
  {
    x2 = x2 + m2;
  }

  x = ( int ) x2;

  return x;
}
//****************************************************************************80

void power_series1 ( int n, double alpha, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_SERIES1 computes a power series for a function G(Z) = (1+F(Z))**ALPHA.
//
//  Discussion:
//
//    The power series for F(Z) is given.
//
//    The form of the power series are:
//
//      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
//
//      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of terms in the power series.
//
//    Input, double ALPHA, the exponent of 1+F(Z) in the definition of G(Z).
//
//    Input, double A[N], the power series coefficients for F(Z).
//
//    Output, double B[N], the power series coefficients for G(Z).
//
{
  int i;
  int j;
  double v;

  for ( j = 1; j <= n; j++ )
  {
    v = 0.0;
    for ( i = 1; i <= j-1; i++ )
    {
      v = v + b[i-1] * a[j-i-1] * ( alpha * ( j - i ) - i );
    }

    b[j-1] = alpha * a[j-1] + v / ( ( double ) j );
  }

  return;
}
//****************************************************************************80

void power_series2 ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_SERIES2 computes the power series for a function G(Z) = EXP(F(Z)) - 1.
//
//  Discussion:
//
//    The power series for F(Z) is given.
//
//    The power series have the form:
//
//      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
//
//      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of terms in the power series.
//
//    Input, double A[N], the power series coefficients for F(Z).
//
//    Output, double B[N], the power series coefficients for G(Z).
//
{
  int i;
  int j;
  double v;

  for ( j = 1; j <= n; j++ )
  {
    v = 0.0;

    for ( i = 1; i <= j-1; i++ )
    {
      v = v + b[i-1] * a[j-i-1] * ( double ) ( j - i );
    }

    b[j-1] = a[j-1] + v / ( double ) j;
  }

  return;
}
//****************************************************************************80

void power_series3 ( int n, double a[], double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_SERIES3 computes the power series for a function H(Z) = G(F(Z)).
//
//  Discussion:
//
//    The power series for F and G are given.
//
//    We assume that
//
//      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
//      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
//      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of terms in the power series.
//
//    Input, double A[N], the power series for F
//
//    Input, double B[N], the power series for G.
//
//    Output, double C[N], the power series for H.
//
{
  int i;
  int iq;
  int j;
  int m;
  double r;
  double v;
  double *work;

  work = new double[n];

  for ( i = 0; i < n; i++ )
  {
    work[i] = b[0] * a[i];
  }
//
//  Search for IQ, the index of the first nonzero entry in A.
//
  iq = 0;

  for ( i = 1; i <= n; i++ )
  {
    if ( a[i-1] != 0.0 )
    {
      iq = i;
      break;
    }
  }

  if ( iq != 0 )
  {
    m = 1;

    for ( ; ; )
    {
      m = m + 1;

      if ( n < m * iq )
      {
        break;
      }

      if ( b[m-1] == 0.0 )
      {
        continue;
      }

      r = b[m-1] * pow ( a[iq-1], m );
      work[m*iq-1] = work[m*iq-1] + r;

      for ( j = 1; j <= n-m*iq; j++ )
      {
        v = 0.0;
        for ( i = 1; i <= j-1; i++ )
        {
          v = v + c[i-1] * a[j-i+iq-1] * ( double ) ( m * ( j - i ) - i );
        }

        c[j-1] = ( ( double ) m * a[j-1] + v / ( double ) j ) / a[iq-1];

      }

      for ( i = 1; i <= n-m*iq; i++ )
      {
        work[i+m*iq-1] = work[i+m*iq-1] + c[i-1] * r;
      }
    }
  }

  for ( i = 0; i < n; i++ )
  {
    c[i] = work[i];
  }

  delete [] work;

  return;
}
//****************************************************************************80

void power_series4 ( int n, double a[], double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_SERIES4 computes the power series for a function H(Z) = G ( 1/F(Z) ).
//
//  Discussion:
//
//    POWER_SERIES4 is given the power series for the functions F and G.
//
//    We assume that
//
//      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
//      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
//      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of terms in the power series.
//
//    Input, double A[N], the power series for F.  For this problem, A(1)
//    may not be 0.0.
//
//    Input, double B(N), the power series for G.
//
//    Output, double C(N), the power series for H.
//
{
  int i;
  int l;
  int m;
  double s;
  double t;
  double *work;

  if ( a[0] == 0.0 )
  {
    cerr << "\n";
    cerr << "POWER_SERIES4 - Fatal error!\n";
    cerr << "  First entry of A is zero.\n";
    exit ( 1 );
  }

  work = new double[n];

  t = 1.0;

  for ( i = 0; i < n; i++ )
  {
    t = t / a[0];
    c[i] = b[i] * t;
    work[i] = a[i] * t;
  }

  for ( m = 2; m <= n; m++ )
  {
    s = -work[m-1];
    for ( i = m; i <= n; i++ )
    {
      for ( l = i; l <= n; l++ )
      {
        c[l-1] = c[l-1] + s * c[l-m];
        work[l-1] = work[l-1] + s * work[l-m];
      }
    }
  }

  delete [] work;

  return;
}
//****************************************************************************80

int prime ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME returns any of the first PRIME_MAX prime numbers.
//
//  Discussion:
//
//    PRIME_MAX is 1600, and the largest prime stored is 13499.
//
//    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 95-98.
//
//  Parameters:
//
//    Input, int N, the index of the desired prime number.
//    In general, is should be true that 0 <= N <= PRIME_MAX.
//    N = -1 returns PRIME_MAX, the index of the largest prime available.
//    N = 0 is legal, returning PRIME = 1.
//
//    Output, int PRIME, the N-th prime.  If N is out of range, PRIME
//    is returned as -1.
//
{
# define PRIME_MAX 1600

  int npvec[PRIME_MAX] = {
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541,
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657,
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657,
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553,
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, 
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, 
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, 
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, 
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, 
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, 
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, 
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, 
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, 
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 };

  if ( n == -1 )
  {
    return PRIME_MAX;
  }
  else if ( n == 0 )
  {
    return 1;
  }
  else if ( n <= PRIME_MAX )
  {
    return npvec[n-1];
  }
  else
  {
    cerr << "\n";
    cerr << "PRIME - Fatal error!\n";
    cerr << "  Unexpected input value of n = " << n << "\n";
    exit ( 1 );
  }

  return 0;
# undef PRIME_MAX
}
//****************************************************************************80

void pythag_triple_next ( int &i, int &j, int &a, int &b, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    PYTHAG_TRIPLE_NEXT computes the next Pythagorean triple.
//
//  Example:
//
//     I       J       A       B       C    A^2+B^2 = C^2
//
//     2       1       3       4       5      25
//     3       2       5      12      13     169
//     4       1      15       8      17     289
//     4       3       7      24      25     625
//     5       2      21      20      29     841
//     5       4       9      40      41    1681
//     6       1      35      12      37    1369
//     6       3      27      36      45    2025
//     6       5      11      60      61    3721
//     7       2      45      28      53    2809
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &I, &J, the generators.
//    On first call, set I = J = 0.  On repeated calls, leave I and J
//    at their output values from the previous call.
//
//    Output, int &A, &B, &C, the next Pythagorean triple.
//    A, B, and C are positive integers which have no common factors,
//    and A**2 + B**2 = C**2.
//
{
//
//  I starts at 2 and increases;
//
//  J starts out at 2 if I is odd, or 1 if I is even, increases by 2,
//    but is always less than I.
//
  
  if ( i == 0 && j == 0 )
  {
    i = 2;
    j = 1;
  }
  else if ( j + 2 < i )
  {
    j = j + 2;
  }
  else
  {
    i = i + 1;
    j = ( i % 2 ) + 1;
  }

  a = i * i - j * j;
  b = 2 * i * j;
  c = i * i + j * j;

  return;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
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
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r4_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

double r8_agm ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    R8_AGM finds the arithmetic-geometric mean of two R8's.
//
//  Discussion:
//
//    The AGM of (A,B) is produced by the following iteration:
//
//      (A,B) -> ( (A+B)/2, SQRT(A*B) ).
//
//    The sequence of successive values of (A,B) quickly converge to the AGM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the numbers whose AGM is desired.  A and B should
//    both be non-negative.
//
//    Output, double R8_AGM, the AGM of the two numbers.
//
{
  double a1;
  double a2;
  double b1;
  double b2;
  double tol;

  if ( a < 0.0 )
  {
    return -1.0;
  }

  if ( b < 0.0 )
  {
    return -1.0;
  }

  if ( a == 0.0 || b == 0.0 )
  {
    return 0.0;
  }

  if ( a == b )
  {
    return a;
  }

  tol = r8_epsilon ( ) * ( a + b + 1.0 );

  a1 = a;
  b1 = b;

  for ( ; ; )
  {
    a2 = ( a1 + b1 ) / 2.0;
    b2 = sqrt ( a1 * b1 );

    if ( fabs ( a2 - b2 ) <= tol )
    {
      return ( ( a2 + b2 ) / 2.0 );
    }

    a1 = a2;
    b1 = b2;
  }
}
//****************************************************************************80

double r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the combinatorial coefficient C(N,K).
//
//  Discussion:
//
//    Real arithmetic is used, and C(N,K) is computed directly, via
//    Gamma functions, rather than recursively.
//
//    C(N,K) is the number of distinct combinations of K objects
//    chosen from a set of N distinct objects.  A combination is
//    like a set, in that order does not matter.
//
//    The formula is:
//
//      C(N,K) = N! / ( (N-K)! * K! )
//
//  Example:
//
//    The number of combinations of 2 things chosen from 5 is 10.
//
//    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
//
//    The actual combinations may be represented as:
//
//      (1,2), (1,3), (1,4), (1,5), (2,3),
//      (2,4), (2,5), (3,4), (3,5), (4,5).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the value of N.
//
//    Input, int K, the value of K.
//
//    Output, double R8_CHOOSE, the value of C(N,K)
//
{
  double arg;
  double fack;
  double facn;
  double facnmk;
  double value;

  if ( n < 0 )
  {
    value = 0.0;
  }
  else if ( k == 0 )
  {
    value = 1.0;
  }
  else if ( k == 1 )
  {
    value = ( double ) n;
  }
  else if ( 1 < k && k < n-1 )
  {
    arg = ( double ) ( n + 1 );
    facn = r8_gamma_log ( arg );

    arg = ( double ) ( k + 1 );
    fack = r8_gamma_log ( arg );

    arg = ( double ) ( n - k + 1 );
    facnmk = r8_gamma_log ( arg );

    value = ( int ) ( 0.5 + exp ( facn - fack - facnmk ) );
  }
  else if ( k == n-1 )
  {
    value = ( double ) n;
  }
  else if ( k == n )
  {
    value = 1.0;
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 round off unit.
//
//  Discussion:
//
//    R8_EPSILON is a number R which is a power of 2 with the property that,
//    to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the round-off unit.
//
{
  double r;

  r = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
//****************************************************************************80

double r8_fall ( double x, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FALL computes the falling factorial function [X]_N.
//
//  Discussion:
//
//    Note that the number of "injections" or 1-to-1 mappings from
//    a set of N elements to a set of M elements is [M]_N.
//
//    The number of permutations of N objects out of M is [M]_N.
//
//    Moreover, the Stirling numbers of the first kind can be used
//    to convert a falling factorial into a polynomial, as follows:
//
//      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
//
//    The formula is:
//
//      [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the falling factorial function.
//
//    Input, int N, the order of the falling factorial function.
//    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
//    negative, a "rising" factorial will be computed.
//
//    Output, double R8_FALL, the value of the falling factorial function.
//
{
  int i;
  double value;

  value = 1.0;

  if ( 0 < n )
  {
    for ( i = 1; i <= n; i++ )
    {
      value = value * x;
      x = x - 1.0;
    }
  }
  else if ( n < 0 )
  {
    for ( i = -1; n <= i; i-- )
    {
      value = value * x;
      x = x + 1.0;
    }
  }

  return value;
}
//****************************************************************************80

double r8_gamma_log ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
//
//  Discussion:
//
//    Computation is based on an algorithm outlined in references 1 and 2.
//    The program uses rational functions that theoretically approximate
//    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
//    approximation for 12 < X is from reference 3, while approximations
//    for X < 12.0 are similar to those in reference 1, but are unpublished.
//
//    The accuracy achieved depends on the arithmetic system, the compiler,
//    intrinsic functions, and proper selection of the machine-dependent
//    constants.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz
//    C++ version by John Burkardt
//
//  Reference:
//
//    William Cody, Kenneth Hillstrom,
//    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
//    Mathematics of Computation,
//    Volume 21, Number 98, April 1967, pages 198-203.
//
//    Kenneth Hillstrom,
//    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
//    May 1969.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maely, Charles Mesztenyi,
//    John Rice, Henry Thatcher, Christop Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the Gamma function.  X must be positive.
//
//    Output, double R8_GAMMA_LOG, the logarithm of the Gamma function of X.
//    If X <= 0.0, or if overflow would occur, the program returns the
//    value XINF, the largest representable double precision number.
//
//
//  Explanation of machine-dependent constants
//
//  BETA   - radix for the real number representation.
//
//  MAXEXP - the smallest positive power of BETA that overflows.
//
//  XBIG   - largest argument for which LN(GAMMA(X)) is representable
//           in the machine, i.e., the solution to the equation
//             LN(GAMMA(XBIG)) = BETA**MAXEXP.
//
//  FRTBIG - Rough estimate of the fourth root of XBIG
//
//
//  Approximate values for some important machines are:
//
//                            BETA      MAXEXP         XBIG
//
//  CRAY-1        (S.P.)        2        8191       9.62E+2461
//  Cyber 180/855
//    under NOS   (S.P.)        2        1070       1.72E+319
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)        2         128       4.08E+36
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)        2        1024       2.55D+305
//  IBM 3033      (D.P.)       16          63       4.29D+73
//  VAX D-Format  (D.P.)        2         127       2.05D+36
//  VAX G-Format  (D.P.)        2        1023       1.28D+305
//
//
//                           FRTBIG
//
//  CRAY-1        (S.P.)   3.13E+615
//  Cyber 180/855
//    under NOS   (S.P.)   6.44E+79
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)   1.42E+9
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)   2.25D+76
//  IBM 3033      (D.P.)   2.56D+18
//  VAX D-Format  (D.P.)   1.20D+9
//  VAX G-Format  (D.P.)   1.89D+76
//
{
  double c[7] = {
    -1.910444077728E-03, 
     8.4171387781295E-04, 
    -5.952379913043012E-04, 
     7.93650793500350248E-04, 
    -2.777777777777681622553E-03, 
     8.333333333333333331554247E-02, 
     5.7083835261E-03 };
  double corr;
  double d1 = - 5.772156649015328605195174E-01;
  double d2 =   4.227843350984671393993777E-01;
  double d4 =   1.791759469228055000094023E+00;
  double frtbig = 1.42E+09;
  int i;
  double p1[8] = {
    4.945235359296727046734888E+00, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = {
    4.974607845568932035012064E+00, 
    5.424138599891070494101986E+02, 
    1.550693864978364947665077E+04, 
    1.847932904445632425417223E+05, 
    1.088204769468828767498470E+06, 
    3.338152967987029735917223E+06, 
    5.106661678927352456275255E+06, 
    3.074109054850539556250927E+06 };
  double p4[8] = {
    1.474502166059939948905062E+04, 
    2.426813369486704502836312E+06, 
    1.214755574045093227939592E+08, 
    2.663432449630976949898078E+09, 
    2.940378956634553899906876E+010,
    1.702665737765398868392998E+011,
    4.926125793377430887588120E+011, 
    5.606251856223951465078242E+011 };
  double pnt68 = 0.6796875E+00;
  double q1[8] = {
    6.748212550303777196073036E+01, 
    1.113332393857199323513008E+03, 
    7.738757056935398733233834E+03, 
    2.763987074403340708898585E+04, 
    5.499310206226157329794414E+04, 
    6.161122180066002127833352E+04, 
    3.635127591501940507276287E+04, 
    8.785536302431013170870835E+03 };
  double q2[8] = {
    1.830328399370592604055942E+02, 
    7.765049321445005871323047E+03, 
    1.331903827966074194402448E+05, 
    1.136705821321969608938755E+06, 
    5.267964117437946917577538E+06, 
    1.346701454311101692290052E+07, 
    1.782736530353274213975932E+07, 
    9.533095591844353613395747E+06 };
  double q4[8] = {
    2.690530175870899333379843E+03, 
    6.393885654300092398984238E+05, 
    4.135599930241388052042842E+07, 
    1.120872109616147941376570E+09, 
    1.488613728678813811542398E+010, 
    1.016803586272438228077304E+011, 
    3.417476345507377132798597E+011, 
    4.463158187419713286462081E+011 };
  double res;
  double sqrtpi = 0.9189385332046727417803297E+00;
  double xbig = 4.08E+36;
  double xden;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double xsq;
//
//  Return immediately if the argument is out of range.
//
  if ( x <= 0.0 || xbig < x )
  {
    return r8_huge ( );
  }

  if ( x <= r8_epsilon ( ) )
  {
    res = -log ( x );
  }
  else if ( x <= 1.5 )
  {
    if ( x < pnt68 )
    {
      corr = - log ( x );
      xm1 = x;
    }
    else
    {
      corr = 0.0;
      xm1 = ( x - 0.5 ) - 0.5;
    }

    if ( x <= 0.5 || pnt68 <= x )
    {
      xden = 1.0;
      xnum = 0.0;

      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm1 + p1[i];
        xden = xden * xm1 + q1[i];
      }

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) );
    }
    else
    {
      xm2 = ( x - 0.5 ) - 0.5;
      xden = 1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm2 + p2[i];
        xden = xden * xm2 + q2[i];
      }

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) );

    }
  }
  else if ( x <= 4.0 )
  {
    xm2 = x - 2.0;
    xden = 1.0;
    xnum = 0.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = xnum * xm2 + p2[i];
      xden = xden * xm2 + q2[i];
    }

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
  }
  else if ( x <= 12.0 )
  {
    xm4 = x - 4.0;
    xden = - 1.0;
    xnum = 0.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = xnum * xm4 + p4[i];
      xden = xden * xm4 + q4[i];
    }

    res = d4 + xm4 * ( xnum / xden );
  }
  else
  {
    res = 0.0;

    if ( x <= frtbig )
    {

      res = c[6];
      xsq = x * x;

      for ( i = 0; i < 6; i++ )
      {
        res = res / xsq + c[i];
      }
    }

    res = res / x;
    corr = log ( x );
    res = res + sqrtpi - 0.5 * corr;
    res = res + x * ( corr - 1.0 );
  }

  return res;
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
//    HUGE_VAL is the largest representable legal real number, and is usually
//    defined in math.h, or sometimes in stdlib.h.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" real value.
//
{
  return HUGE_VAL;
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the integer that is nearest to a double value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the real number.
//
//    Output, int R8_NINT, the nearest integer.
//
{
  double d;
  int i;

  i = int ( x );
  d = x - i;

  if ( fabs ( d ) <= 0.5 )
  {
    return i;
  }
  else if ( x < i ) 
  {
    return (i-1);
  }
  else
  {
    return (i+1);
  }
}
//****************************************************************************80

double r8_pi ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PI returns the value of PI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_PI, the value of PI.
//
{
  return ( double ) 3.14159265358979323846264338327950288419716939937510;
}
//****************************************************************************80

double r8_rise ( double x, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_RISE computes the rising factorial function [X]^N.
//
//  Discussion:
//
//    [X}^N = X * ( X + 1 ) * ( X + 2 ) * ... * ( X + N - 1 ).
//
//    Note that the number of ways of arranging N objects in M ordered
//    boxes is [M}^N.  (Here, the ordering in each box matters).  Thus,
//    2 objects in 2 boxes have the following 6 possible arrangements:
//
//      -/12, 1/2, 12/-, -/21, 2/1, 21/-.
//
//    Moreover, the number of non-decreasing maps from a set of
//    N to a set of M ordered elements is [M]^N / N!.  Thus the set of
//    nondecreasing maps from (1,2,3) to (a,b,c,d) is the 20 elements:
//
//      aaa, abb, acc, add, aab, abc, acd, aac, abd, aad
//      bbb, bcc, bdd, bbc, bcd, bbd, ccc, cdd, ccd, ddd.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the rising factorial function.
//
//    Input, int N, the order of the rising factorial function.
//    If N = 0, RISE = 1, if N = 1, RISE = X.  Note that if N is
//    negative, a "falling" factorial will be computed.
//
//    Output, double R8_RISE, the value of the rising factorial function.
//
{
  int i;
  double value;

  value = 1.0;

  if ( 0 < n )
  {
    for ( i = 1; i <= n; i++ )
    {
      value = value * x;
      x = x + 1.0;
    }
  }
  else if ( n < 0 )
  {
    for ( i = -1; n <= i; i-- )
    {
      value = value * x;
      x = x - 1.0;
    }
  }

  return value;
}
//****************************************************************************80

void r8_swap ( double &x, double &y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SWAP switches two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, double &X, &Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  double z;

  z = x;
  x = y;
  y = z;
 
  return;
}
//****************************************************************************80

void r8_to_cfrac ( double r, int n, int a[], int p[], int q[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_CFRAC converts a double value to a continued fraction.
//
//  Discussion:
//
//    The routine is given a double number R.  It computes a sequence of
//    continued fraction approximations to R, returning the results as
//    simple fractions of the form P(I) / Q(I).
//
//  Example:
//
//    X = 2 * PI
//    N = 7
//
//    A = [ *, 6,  3,  1,  1,   7,   2,    146,      3 ]
//    P = [ 1, 6, 19, 25, 44, 333, 710, 103993, 312689 ]
//    Q = [ 0, 1,  3,  4,  7,  53, 113,  16551,  49766 ]
//
//    (This ignores roundoff error, which will cause later terms to differ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Norman Richert,
//    Strang's Strange Figures,
//    American Mathematical Monthly,
//    Volume 99, Number 2, February 1992, pages 101-107.
//
//  Parameters:
//
//    Input, double R, the double value.
//
//    Input, int N, the number of convergents to compute.
//
//    Output, int A[N+1], the partial quotients.
//
//    Output, int P[N+2], Q[N+2], the numerators and denominators
//    of the continued fraction approximations.
//
{
  int i;
  double r_copy;
  double *x;

  if ( r == 0.0 )
  {
    for ( i = 0; i <= n; i++ )
    {
      a[i] = 0;
    }
    for ( i = 0; i <= n+1; i++ )
    {
      p[i] = 0;
    }
    for ( i = 0; i <= n+1; i++ )
    {
      q[i] = 0;
    }
    return;
  }

  x = new double[n+1];

  r_copy = fabs ( r );

  p[0] = 1;
  q[0] = 0;

  p[1] = ( int ) r_copy;
  q[1] = 1;
  x[0] = r_copy;
  a[0] = ( int ) x[0];

  for ( i = 1; i <= n; i++ )
  {
    x[i] = 1.0 / ( x[i-1] - ( double ) a[i-1] );
    a[i] = ( int ) x[i];
    p[i+1] = a[i] * p[i] + p[i-1];
    q[i+1] = a[i] * q[i] + q[i-1];
  }

  if ( r < 0.0 )
  {
    for ( i = 0; i <= n+1; i++ )
    {
      p[i] = -p[i];
    }
  }

  delete [] x;

  return;
}
//****************************************************************************80

void r8_to_dec ( double dval, int dec_digit, int &mantissa, int &exponent )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_DEC converts a double quantity to a decimal representation.
//
//  Discussion:
//
//    Given the double value DVAL, the routine computes integers
//    MANTISSA and EXPONENT so that it is approximatelytruethat:
//
//      DVAL = MANTISSA * 10 ** EXPONENT
//
//    In particular, only DEC_DIGIT digits of DVAL are used in constructing the
//    representation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double DVAL, the value whose decimal representation
//    is desired.
//
//    Input, int DEC_DIGIT, the number of decimal digits to use.
//
//    Output, int &MANTISSA, &EXPONENT, the approximate decimal 
//    representation of DVAL.
//
{
  double mantissa_double;
  double ten1;
  double ten2;
//
//  Special cases.
//
  if ( dval == 0.0 )
  {
    mantissa = 0;
    exponent = 0;
    return;
  }
//
//  Factor DVAL = MANTISSA_DOUBLE * 10**EXPONENT
//
  mantissa_double = dval;
  exponent = 0;
//
//  Now normalize so that 
//  10**(DEC_DIGIT-1) <= ABS(MANTISSA_DOUBLE) < 10**(DEC_DIGIT)
//
  ten1 = pow ( 10.0, dec_digit - 1 );
  ten2 = 10.0 * ten1;

  while ( fabs ( mantissa_double ) < ten1 )
  {
    mantissa_double = mantissa_double * 10.0;
    exponent = exponent - 1;
  }

  while ( ten2 <= fabs ( mantissa_double ) )
  {
    mantissa_double = mantissa_double / 10.0;
    exponent = exponent + 1;
  }
//
//  MANTISSA is the integer part of MANTISSA_DOUBLE, rounded.
//
  mantissa = r8_nint ( mantissa_double );
//
//  Now divide out any factors of ten from MANTISSA.
//
  if ( mantissa != 0 )
  {
    while ( 10 * ( mantissa / 10 ) == mantissa )
    {
      mantissa = mantissa / 10;
      exponent = exponent + 1;
    }
  }

  return;
}
//****************************************************************************80

void r8_to_rat ( double a, int ndig, int &iatop, int &iabot )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_RAT converts a real value to a rational value.
//
//  Discussion:
//
//    The rational value (IATOP/IABOT) is essentially computed by truncating
//    the decimal representation of the real value after a given number of
//    decimal digits.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the real value to be converted.
//
//    Input, int NDIG, the number of decimal digits used.
//
//    Output, int &IATOP, &IABOT, the numerator and denominator
//    of the rational value that approximates A.
//
{
  double factor;
  int i;
  int ifac;
  int itemp;
  int jfac;

  factor = pow ( 10.0, ndig );

  if ( 0 < ndig )
  {
    iabot = 1;
    for ( i = 1; i <= ndig; i++ )
    {
      iabot = iabot * 10;
    }
    iatop = 1;
  }
  else
  {
    iabot = 1;
    iatop = 1;
    for ( i = 1; i <= -ndig; i++ )
    {
      iatop = iatop * 10;
    }
  }

  iatop = r8_nint ( a * factor ) * iatop;
  iabot = iabot;
//
//  Factor out the greatest common factor.
//
  itemp = i4_gcd ( iatop, iabot );

  iatop = iatop / itemp;
  iabot = iabot / itemp;

  return;
}
//****************************************************************************80

double r8_uniform ( double rlo, double rhi, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM returns a random real in a given range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RLO, RHI, the minimum and maximum values.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8_UNIFORM, the randomly chosen value.
//
{
  double t;

  t = r8_uniform_01 ( seed );

  return ( 1.0 - t ) * rlo + t * rhi;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a double precision pseudorandom number.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double r8mat_det ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET finds the determinant of a real N by N matrix.
//
//  Discussion:
//
//    A brute force calculation is made.
//
//    This routine should only be used for small matrices, since this
//    calculation requires the summation of N! products of N numbers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of A.
//
//    Input, double A[N*N], the matrix whose determinant is desired.
//
//    Output, double R8MAT_DET, the determinant of the matrix.
//
{
  double det;
  bool even;
  int i;
  bool more;
  int *perm;
  double term;

  more = false;
  det = 0.0;
  perm = new int[n];

  for ( ; ; )
  {
    perm_next ( n, perm, more, even );

    if ( even )
    {
      term = 1.0;
    }
    else
    {
      term = -1.0;
    }

    for ( i = 1; i <= n; i++ )
    {
      term = term * a[i-1+(perm[i-1]-1)*n];
    }

    det = det + term;

    if ( !more )
    {
      break;
    }
  }

  delete [] perm;

  return det;
}
//****************************************************************************80

void r8mat_perm ( int n, double a[], int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PERM permutes the rows and columns of a square R8MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double A[N*N].
//    On input, the matrix to be permuted.
//    On output, the permuted matrix.
//
//    Input, int P[N], a permutation to be applied to the rows
//    and columns.  P[I] is the new number of row and column I.
//
{
  int i;
  int i1;
  int iopt = 1;
  int isgn;
  double it;
  int j;
  int j1;
  int j2;
  int k;
  int lc;
  int ncycle;
  double temp;

  perm_cycle ( n, p, isgn, ncycle, iopt );

  for ( i = 1; i <= n; i++ )
  {
    i1 = - p[i-1];

    if ( 0 < i1 )
    {

      lc = 0;

      for ( ; ; )
      {
        i1 = p[i1-1];
        lc = lc + 1;

        if ( i1 <= 0 )
        {
          break;
        }

      }

      i1 = i;

      for ( j = 1; j <= n; j++ )
      {
        if ( p[j-1] <= 0 )
        {
          j2 = j;
          k = lc;

          for ( ; ; )
          {
            j1 = j2;
            it = a[i1-1+(j1-1)*n];

            for ( ; ; )
            {
              i1 = abs ( p[i1-1] );
              j1 = abs ( p[j1-1] );

              temp = a[i1-1+(j1-1)*n];
              a[i1-1+(j1-1)*n] = it;
              it = temp;

              if ( j1 != j2 )
              {
                continue;
              }

              k = k - 1;

              if ( i1 == i )
              {
                break;
              }
            }

            j2 = abs ( p[j2-1] );

            if ( k == 0 ) 
            {
              break;
            }
          }
        }
      }
    }
  }
//
//  Restore the positive signs of the data.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = abs ( p[i] );
  }

  return;
}
//****************************************************************************80

void r8mat_perm2 ( int m, int n, double a[], int p[], int q[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PERM2 permutes rows and columns of a rectangular R8MAT, in place.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int M, number of rows in the matrix.
//
//    Input, int N, number of columns in the matrix.
//
//    Input/output, double A[M*N].
//    On input, the matrix to be permuted.
//    On output, the permuted matrix.
//
//    Input, int P[M], the row permutation.  P(I) is the new number of row I.
//
//    Input, int Q[N], the column permutation.  Q(I) is the new number of
//    column I.  
//
{
  int i;
  int i1;
  int is;
  int j;
  int j1;
  int j2;
  int k;
  int lc;
  int nc;
  double t;
  double temp;

  perm_cycle ( m, p, is, nc, 1 );

  perm_cycle ( n, q, is, nc, 1 );

  for ( i = 1; i <= m; i++ )
  {
    i1 = - p[i-1];

    if ( 0 < i1 )
    {
      lc = 0;

      for ( ; ; )
      {
        i1 = p[i1-1];
        lc = lc + 1;

        if ( i1 <= 0 )
        {
          break;
        }

      }

      i1 = i;

      for ( j = 1; j <= n; j++ )
      {
        if ( q[j-1] <= 0 )
        {
          j2 = j;
          k = lc;

          for ( ; ; )
          {
            j1 = j2;
            t = a[i1-1+(j1-1)*n];

            for ( ; ; )
            {
              i1 = abs ( p[i1-1] );
              j1 = abs ( q[j1-1] );

              temp = a[i1-1+(j1-1)*n];
              a[i1-1+(j1-1)*n] = t;
              t = temp;
 

              if ( j1 != j2 )
              {
                continue;
              }

              k = k - 1;

              if ( i1 == i )
              {
                break;
              }
            }

            j2 = abs ( q[j2-1] );

            if ( k == 0 )
            {
              break;
            }
          }
        }
      }
    }
  }

  for ( i = 0; i < m; i++ )
  {
    p[i] = abs ( p[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    q[i] = abs ( q[i] );
  }

  return;
}
//****************************************************************************80

double r8mat_permanent ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PERMANENT computes the permanent of an R8MAT.
//
//  Discussion:
//
//    The permanent function is similar to the determinant.  Recall that
//    the determinant of a matrix may be defined as the sum of all the
//    products:
//
//      S * A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
//
//    where I is any permutation of the columns of the matrix, and S is the
//    sign of the permutation.  By contrast, the permanent function is
//    the (unsigned) sum of all the products
//
//      A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
//
//    where I is any permutation of the columns of the matrix.  The only
//    difference is that there is no permutation sign multiplying each summand.
//
//    Symbolically, then, the determinant of a 2 by 2 matrix
//
//      a b
//      c d
//
//    is a*d-b*c, whereas the permanent of the same matrix is a*d+b*c.
//
//
//    The permanent is invariant under row and column permutations.
//    If a row or column of the matrix is multiplied by S, then the
//      permanent is likewise multiplied by S.
//    If the matrix is square, then the permanent is unaffected by
//      transposing the matrix.
//    Unlike the determinant, however, the permanent does change if
//      one row is added to another, and it is not true that the
//      permanent of the product is the product of the permanents.
//
//
//    Note that if A is a matrix of all 1's and 0's, then the permanent
//    of A counts exactly which permutations hit exactly 1's in the matrix.
//    This fact can be exploited for various combinatorial purposes.
//
//    For instance, setting the diagonal of A to 0, and the offdiagonals
//    to 1, the permanent of A counts the number of derangements of N
//    objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, number of rows and columns in matrix.
//
//    Input, double A[N*N], the matrix whose permanent is desired.
//
//    Output, double R8MAT_PERMANENT, the value of the permanent of A.
//
{
  int i;
  int iadd;
  int *iwork;
  int j;
  bool more;
  int ncard;
  double p;
  double perm;
  double prod;
  double sgn;
  double *work;
  double z;

  more = false;

  iwork = new int[n];
  work = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    work[i-1] = a[i-1+(n-1)*n];
    for ( j = 1; j <= n; j++ )
    {
      work[i-1] = work[i-1] - 0.5 * a[i-1+(j-1)*n];
    }
  }

  p = 0.0;
  sgn = -1.0;

  for ( ; ; )
  {
    sgn = -sgn;

    subset_gray_next ( n-1, iwork, more, ncard, iadd );

    if ( ncard != 0 )
    {
      z = ( double ) ( 2 * iwork[iadd-1] - 1 );
      for ( i = 1; i <= n; i++ )
      {
        work[i-1] = work[i-1] + z * a[i-1+(iadd-1)*n];
      }
    }

    prod = 1.0;
    for ( i = 0; i < n; i++ )
    {
      prod = prod * work[i];
    }
    p = p + sgn * prod;

    if ( !more )
    {
      break;
    }
  }

  delete [] iwork;
  delete [] work;

  perm = p * ( double ) ( 4 * ( n % 2 ) - 2 );

  return perm;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT, with an optional title.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
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
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
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

void r8poly ( int n, double a[], double x0, int iopt, double &val )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY performs operations on R8POLY's in power or factorial form.
//
//  Discussion:
//
//    The power sum form of a polynomial is
//
//      P(X) = A1 + A2*X + A3*X**2 + ... + (AN+1)*X**N
//
//    The Taylor expansion at C has the form
//
//      P(X) = A1 + A2*(X-C) + A3*(X-C)**2+... + (AN+1)*(X-C)**N
//
//    The factorial form of a polynomial is
//
//      P(X) = A1 + A2*X + A3*(X)*(X-1) + A4*(X)*(X-1)*(X-2) + ...
//        + (AN+1)*(X)*(X-1)*...*(X-N+1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of coefficients in the polynomial
//    (in other words, the polynomial degree + 1)
//
//    Input/output, double A[N], the coefficients of the polynomial.  Depending
//    on the option chosen, these coefficients may be overwritten by those
//    of a different form of the polynomial.
//
//    Input, double X0, for IOPT = -1, 0, or positive, the value of the
//    argument at which the polynomial is to be evaluated, or the
//    Taylor expansion is to be carried out.
//
//    Input, int IOPT, a flag describing which algorithm is to
//    be carried out:
//    -3: Reverse Stirling.  Input the coefficients of the polynomial in
//    factorial form, output them in power sum form.
//    -2: Stirling.  Input the coefficients in power sum
//    form, output them in factorial form.
//    -1: Evaluate a polynomial which has been input
//    in factorial form.
//    0:  Evaluate a polynomial input in power sum form.
//    1 or more:  Given the coefficients of a polynomial in
//    power sum form, compute the first IOPT coefficients of
//    the polynomial in Taylor expansion form.
//
//    Output, double &VAL, for IOPT = -1 or 0, the value of the
//    polynomial at the point X0.
//
{
  double eps;
  int i;
  int m;
  int n1;
  double w;
  double z;

  n1 = i4_min ( n, iopt );
  n1 = i4_max ( 1, n1 );

  if ( iopt < -1 )
  {
    n1 = n;
  }

  eps = ( double ) ( i4_max ( -iopt, 0 ) % 2 );

  w = - ( double ) n * eps;

  if ( -2 < iopt )
  {
    w = w + x0;
  }

  for ( m = 1; m <= n1; m++ )
  {
    val = 0.0;
    z = w;

    for ( i = m; i <= n; i++ )
    {
      z = z + eps;
      val = a[n+m-i-1] + z * val;
      if ( iopt != 0 && iopt != -1 )
      {
        a[n+m-i-1] = val;
      }
    }

    if ( iopt < 0 )
    {
      w = w + 1.0;
    }
  }

  return;
}
//****************************************************************************80

int r8poly_degree ( int na, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_DEGREE returns the degree of an R8POLY.
//
//  Discussion:
//
//    The degree of a polynomial is the index of the highest power
//    of X with a nonzero coefficient.
//
//    The degree of a constant polynomial is 0.  The degree of the
//    zero polynomial is debatable, but this routine returns the
//    degree as 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the dimension of A.
//
//    Input, double A[NA+1], the coefficients of the polynomials.
//
//    Output, int R8POLY_DEGREE, the degree of the polynomial.
//
{
  int degree;

  degree = na;

  while ( 0 < degree )
  {
    if ( a[degree] != 0.0 )
    {
      return degree;
    }
    degree = degree - 1;
  }
  return degree;
}
//****************************************************************************80

void r8poly_div ( int na, double a[], int nb, double b[], int &nq, double q[], 
  int &nr, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_DIV computes the quotient and remainder of two R8POLY's.
//
//  Discussion:
//
//    The polynomials are assumed to be stored in power sum form.
//
//    The power sum form is:
//
//      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the dimension of A.
//
//    Input, double A[NA+1], the coefficients of the polynomial to be divided.
//
//    Input, int NB, the dimension of B.
//
//    Input, double B[NB+1], the coefficients of the divisor polynomial.
//
//    Output, int &NQ, the degree of Q.
//    If the divisor polynomial is zero, NQ is returned as -1.
//
//    Output, double Q[NA+1-NB], contains the quotient of A/B.
//    If A and B have full degree, Q should be dimensioned Q(0:NA-NB).
//    In any case, Q(0:NA) should be enough.
//
//    Output, int &NR, the degree of R.
//    If the divisor polynomial is zero, NR is returned as -1.
//
//    Output, double R[NB], contains the remainder of A/B.
//    If B has full degree, R should be dimensioned R(0:NB-1).
//    Otherwise, R will actually require less space.
//
{
  double *a2;
  int i;
  int j;
  int na2;
  int nb2;

  na2 = r8poly_degree ( na, a );
  nb2 = r8poly_degree ( nb, b );

  if ( b[nb2] == 0.0 )
  {
    nq = -1;
    nr = -1;
    return;
  }

  a2 = new double[na+1];
  for ( i = 0; i <= na; i++ )
  {
    a2[i] = a[i];
  }

  nq = na2 - nb2;
  nr = nb2 - 1;

  for ( i = nq; 0 <= i; i-- )
  {
    q[i] = a2[i+nb2] / b[nb2];
    a2[i+nb2] = 0.0;
    for ( j = 0; j <= nb2-1; j++ )
    {
      a2[i+j] = a2[i+j] - q[i] * b[j];
    }
  }

  for ( i = 0; i <= nr; i++ )
  {
    r[i] = a2[i];
  }

  delete [] a2;

  return;
}
//****************************************************************************80

void r8poly_f2p ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_F2P converts a double polynomial from factorial form to power sum form.
//
//  Discussion:
//
//    The (falling) factorial form is
//
//      p(x) =   a(1)
//             + a(2) * x
//             + a(3) * x*(x-1)
//             ...
//             + a(n) * x*(x-1)*...*(x-(n-2))
//
//    The power sum form is
//
//      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input/output, double A[N], on input, the polynomial
//    coefficients in factorial form.  On output, the polynomial
//    coefficients in power sum form.
//
{
  int i;
  int m;
  double val;
  double w;
  double z;

  w = - ( double ) n;

  for ( m = 1; m <= n; m++ )
  {
    val = 0.0;
    z = w;

    for ( i = m; i <= n; i++ )
    {
      z = z + 1.0;
      val = a[n+m-i-1] + z * val;
      a[n+m-i-1] = val;
    }
    w = w + 1.0;
  }
  return;
}
//****************************************************************************80

double r8poly_fval ( int n, double a[], double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_FVAL evaluates a double polynomial in factorial form.
//
//  Discussion:
//
//    The (falling) factorial form of a polynomial is:
//
//      p(x) = a(1)
//           + a(2)  *x
//           + a(3)  *x*(x-1)
//           +...
//           + a(n-1)*x*(x-1)*(x-2)...*(x-(n-3))
//           + a(n)  *x*(x-1)*(x-2)...*(x-(n-3))*(x-(n-2))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input, double A[N], the coefficients of the polynomial.
//    A(1) is the constant term.
//
//    Input, double X, the point at which the polynomial is to be evaluated.
//
//    Output, double R8POLY_FVAL, the value of the polynomial at X.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = 1; i <= n; i++ )
  {
    value = a[n-i] + ( x - n + i ) * value;
  }

  return value;
}
//****************************************************************************80

void r8poly_mul ( int na, double a[], int nb, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_MUL computes the product of two double polynomials A and B.
//
//  Discussion:
//
//    The polynomials are in power sum form.
//
//    The power sum form is:
//
//      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the dimension of A.
//
//    Input, double A[NA+1], the coefficients of the first polynomial factor.
//
//    Input, int NB, the dimension of B.
//
//    Input, double B[NB+1], the coefficients of the second polynomial factor.
//
//    Output, double C[NA+NB+1], the coefficients of A * B.
//
{
  double *d;
  int i;
  int j;

  d = new double[na+nb+1];

  for ( i = 0; i <= na + nb; i++ )
  {
    d[i] = 0.0;
  }

  for ( i = 0; i <= na; i++ )
  {
    for ( j = 0; j <= nb; j++ )
    {
      d[i+j] = d[i+j] + a[i] * b[j];
    }
  }

  for ( i = 0; i <= na + nb; i++ )
  {
    c[i] = d[i];
  }

  delete [] d;

  return;
}
//****************************************************************************80

void r8poly_n2p ( int n, double a[], double xarray[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_N2P converts a double polynomial from Newton form to power sum form.
//
//  Discussion:
//
//    This is done by shifting all the Newton abscissas to zero.
//
//    Actually, what happens is that the abscissas of the Newton form
//    are all shifted to zero, which means that A is the power sum
//    polynomial description and A, XARRAY is the Newton polynomial
//    description.  It is only because all the abscissas are shifted to
//    zero that A can be used as both a power sum and Newton polynomial
//    coefficient array.
//
//    The Newton form of a polynomial is described by an array of N coefficients
//    A and N abscissas X:
//
//      p(x) =   a(1)
//             + a(2) * (x-x(1))
//             + a(3) * (x-x(1)) * (x-x(2))
//             ...
//             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
//
//    X(N) does not occur explicitly in the formula for the evaluation of p(x),
//    although it is used in deriving the coefficients A.
//
//    The power sum form of a polynomial is:
//
//      p(x) = a(1) + a(2)*x + ... + a(n-1)*x**(n-2) + a(n)*x**(n-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input/output, double A[N].  On input, the coefficients
//    of the polynomial in Newton form, and on output, the coefficients
//    in power sum form.
//
//    Input/output, double XARRAY[N].  On input, the abscissas of
//    the Newton form of the polynomial.  On output, these values
//    have all been set to zero.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    r8poly_nx ( n, a, xarray, 0.0 );
  }

  return;
}
//****************************************************************************80

double r8poly_nval ( int n, double a[], double xarray[], double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_NVAL evaluates a double polynomial in Newton form.
//
//  Definition:
//
//    The Newton form of a polynomial is;
//
//      p(x) = a(1)
//           + a(2)  *(x-x1)
//           + a(3)  *(x-x1)*(x-x2)
//           +...
//           + a(n-1)*(x-x1)*(x-x2)*(x-x3)...*(x-x(n-2))
//           + a(n)  *(x-x1)*(x-x2)*(x-x3)...*(x-x(n-2))*(x-x(n-1))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input, double A[N], the coefficients of the polynomial.
//    A(1) is the constant term.
//
//    Input, double XARRAY[N-1], the N-1 points X which are part
//    of the definition of the polynomial.
//
//    Input, double X, the point at which the polynomial is to be evaluated.
//
//    Output, double R8POLY_NVAL, the value of the polynomial at X.
//
{
  int i;
  double value;

  value = a[n-1];

  for ( i = n-2; 0 <= i; i-- )
  {
    value = a[i] + ( x - xarray[i] ) * value;
  }

  return value;
}
//****************************************************************************80

void r8poly_nx ( int n, double a[], double xarray[], double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_NX replaces one of the base points in a polynomial in Newton form.
//
//  Discussion:
//
//    The Newton form of a polynomial is described by an array of N coefficients
//    A and N abscissas X:
//
//      p(x) =   a(1)
//             + a(2) * (x-x(1))
//             + a(3) * (x-x(1)) * (x-x(2))
//             ...
//             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
//
//    X(N) does not occur explicitly in the formula for the evaluation of p(x),
//    although it is used in deriving the coefficients A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//   Input, int N, the dimension of A.
//
//   Input/output, double A[N], the polynomial coefficients of the Newton form.
//
//   Input/output, double XARRAY[N], the set of abscissas that
//   are part of the Newton form of the polynomial.  On output,
//   the abscissas have been shifted up one index, so that
//   the first location now holds X, and the original value
//   of the last entry is discarded.
//
//   Input, double X, the new point to be shifted into XARRAY.
//
{
  int i;

  for ( i = n - 2; 0 <= i; i-- )
  {
    a[i] = a[i] + ( x - xarray[i] ) * a[i+1];
  }

  for ( i = n - 1; 0 < i; i-- )
  {
    xarray[i] = xarray[i-1];
  }

  xarray[0] = x;

  return;
}
//****************************************************************************80

void r8poly_p2f ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_P2F converts a double polynomial from power sum form to factorial form.
//
//  Discussion:
//
//    The power sum form is
//
//      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
//
//    The (falling) factorial form is
//
//      p(x) =   a(1)
//             + a(2) * x
//             + a(3) * x*(x-1)
//             ...
//             + a(n) * x*(x-1)*...*(x-(n-2))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input/output, double A[N], on input, the polynomial
//    coefficients in the power sum form, on output, the polynomial
//    coefficients in factorial form.
//
{
  int i;
  int m;
  double val;

  for ( m = 1; m <= n; m++ )
  {
    val = 0.0;
    for ( i = m; i <= n; i++ )
    {
      val = a[n+m-i-1] + ( double ) ( m - 1 ) * val;
      a[n+m-i-1] = val;
    }
  }

  return;
}
//****************************************************************************80

void r8poly_p2n ( int n, double a[], double xarray[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_P2N converts a polynomial from power sum form to Newton form.
//
//  Discussion:
//
//    This is done by shifting all the Newton abscissas from zero.
//
//    The power sum form of a polynomial is:
//
//      p(x) = a(1) + a(2) * x + ... + a(n-1) * x**(n-2) + a(n) * x**(n-1)
//
//    The Newton form of a polynomial is described by an array of N coefficients
//    A and N abscissas X:
//
//      p(x) =   a(1)
//             + a(2) * (x-x(1))
//             + a(3) * (x-x(1)) * (x-x(2))
//             ...
//             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
//
//    X(N) does not occur explicitly in the formula for the evaluation of p(x),
//    although it is used in deriving the coefficients A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input/output, double A[N].  On input, the coefficients
//    of the polynomial in power sum form, and on output, the
//    coefficients in Newton form.
//
//    Input, double XARRAY[N].  On input, the desired abscissas of
//    the Newton form of the polynomial.
//
{
  int i;
  double *work;
  double x;

  work = new double[n];

  for ( i = 0; i < n; i++ )
  {
    work[i] = 0.0;
  }

  for ( i = n - 1; 0 <= i; i-- )
  {
    x = xarray[i];
    r8poly_nx ( n, a, work, x );
  }

  delete [] work;

  return;
}
//****************************************************************************80

void r8poly_p2t ( int n, double a[], double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_P2T converts a real polynomial from power sum form to Taylor form.
//
//  Discussion:
//
//    The power sum form is
//
//      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
//
//    The Taylor form is
//
//      p(x) =   a(1)
//             + a(2) * (x-x0)
//             + a(3) * (x-x0)**2
//             ...
//             + a(n) * (x-x0)**(n-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input/output, double A[N], on input, the coefficients in
//    power sum form, and on output, the coefficients in Taylor form.
//
//    Input, double X, the point at which the Taylor form of the
//    polynomial is to be based.
//
{
  int i;
  int m;
  double val;

  for ( m = 1; m <= n; m++ )
  {
    val = 0.0;
    for ( i = m; i <= n; i++ )
    {
      val = a[n+m-i-1] + x * val;
      a[n+m-i-1] = val;
    }
  }

  return;
}
//****************************************************************************80

void r8poly_power ( int na, double a[], int p, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_POWER computes a positive integer power of a polynomial.
//
//  Discussion:
//
//    The power sum form is:
//
//      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the dimension of A.
//
//    Input, double A[NA+1], the polynomial to be raised to the power.
//
//    Input, int P, the nonnegative power to which A is raised.
//
//    Output, double B[P*NA+1], the power of the polynomial.
//
{
  int i;
  int j;
  int nonzer;
//
//  Zero out B.
//
  for ( i = 0; i <= p * na; i++ )
  {
    b[i] = 0.0;
  }
//
//  Search for the first nonzero element in A.
//
  nonzer = 0;

  for ( i = 0; i <= na; i++ )
  {
    if ( a[i] != 0.0 )
    {
      nonzer = i;
      break;
    }
  }

  if ( nonzer == 0 )
  {
    return;
  }

  b[0] = pow ( a[nonzer], p );

  for ( i = 1; i <= p * ( na - nonzer ); i++ )
  {
    if ( i + nonzer <= na )
    {
      b[i] = ( double ) ( i * p ) * b[0] * a[i+nonzer];
    }
    else
    {
      b[i] = 0.0;
    }

    for ( j = 1; j <= i-1; j++ )
    {
      if ( j+nonzer <= na )
      {
        b[i] = b[i] - ( double ) ( i - j ) * a[j+nonzer] * b[i-j];
      }

      if ( i-j+nonzer <= na )
      {
        b[i] = b[i] + ( double ) ( i - j ) * ( double ) p * b[j] * a[i-j+nonzer];
      }
    }

    b[i] = b[i] / ( ( double ) i  * a[nonzer] );
  }
//
//  Shift B up.
//
  for ( i = p*nonzer; i <= p*na; i++ )
  {
    b[i] = b[i-p*nonzer];
  }

  for ( i = 0; i <= p * nonzer-1; i++ )
  {
    b[i] = 0.0;
  }

  return;
}
//****************************************************************************80

void r8poly_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_PRINT prints out a polynomial.
//
//  Discussion:
//
//    The power sum form is:
//
//      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input, double A[N+1], the polynomial coefficients.
//    A(0) is the constant term and
//    A(N) is the coefficient of X**N.
//
//    Input, string TITLE, a title.
//
{
  int i;
  double mag;
  int n2;
  char plus_minus;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  n2 = r8poly_degree ( n, a );

  if ( a[n2] < 0.0 )
  {
    plus_minus = '-';
  }
  else
  {
    plus_minus = ' ';
  }

  mag = fabs ( a[n2] );

  if ( 2 <= n2 )
  {
    cout << "  p(x) = " << plus_minus << mag << " * x^" << n2 << "\n";
  }
  else if ( n2 == 1 )
  {
    cout << "  p(x) = " << plus_minus << mag << " * x" << "\n";
  }
  else if ( n2 == 0 )
  {
    cout << "  p(x) = " << plus_minus << mag << "\n";
  }

  for ( i = n2-1; 0 <= i; i-- )
  {
    if ( a[i] < 0.0 )
    {
      plus_minus = '-';
    }
    else
    {
      plus_minus = '+';
    }

    mag = fabs ( a[i] );

    if ( mag != 0.0 )
    {
      if ( 2 <= i )
      {
        cout << "         " << plus_minus << mag << " * x^" << i << "\n";
      }
      else if ( i == 1 )
      {
        cout << "         " << plus_minus << mag << " * x" << "\n";
      }
      else if ( i == 0 )
      {
        cout << "         " << plus_minus << mag << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

double r8poly_pval ( int n, double a[], double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_PVAL evaluates a real polynomial in power sum form.
//
//  Discussion:
//
//    The power sum form is:
//
//      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input, double A[N+1], the coefficients of the polynomial.
//    A(0) is the constant term.
//
//    Input, double X, the point at which the polynomial is to be evaluated.
//
//    Output, double R8POLY_VAL, the value of the polynomial at X.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = n; 0 <= i; i-- )
  {
    value = value * x + a[i];
  }

  return value;
}
//****************************************************************************80

void r8poly_t2p ( int n, double a[], double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_T2P converts a real polynomial from Taylor form to power sum form
//
//  Discussion:
//
//    The Taylor form is
//
//      p(x) =   a(1)
//             + a(2) * (x-x0)
//             + a(3) * (x-x0)**2
//             ...
//             + a(n) * (x-x0)**(n-1)
//
//    The power sum form is
//
//      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input/output, double A[N].  On input, the coefficients in Taylor form,
//    and on output, the coefficients in power sum form.
//
//    Input, double X, the point at which the Taylor form polynomial is based.
//
{
  int i;
  int j;

  for ( i = n; 1 <= i; i-- )
  {
    for ( j = i; j <= n-1; j++ )
    {
      a[j-1] = a[j-1] - a[j] * x;
    }
  }

  return;
}
//****************************************************************************80

void r8vec_backtrack ( int n, int maxstack, int stack[], double x[], int &indx, 
  int &k, int &nstack, int ncan[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BACKTRACK supervises a backtrack search for a real vector.
//
//  Discussion:
//
//    The routine tries to construct a real vector one index at a time,
//    using possible candidates as supplied by the user.
//
//    At any time, the partially constructed vector may be discovered to be
//    unsatisfactory, but the routine records information about where the
//    last arbitrary choice was made, so that the search can be
//    carried out efficiently, rather than starting out all over again.
//
//    First, call the routine with INDX = 0 so it can initialize itself.
//
//    Now, on each return from the routine, if INDX is:
//      1, you've just been handed a complete candidate vector;
//         Admire it, analyze it, do what you like.
//      2, please determine suitable candidates for position X(K).
//         Return the number of candidates in NCAN(K), adding each
//         candidate to the end of STACK, and increasing NSTACK.
//      3, you're done.  Stop calling the routine;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 July 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of positions to be filled in the vector.
//
//    Input, int MAXSTACK, the maximum length of the stack.
//
//    Input, double STACK[MAXSTACK], a list of all current candidates for
//    all positions 1 through K.
//
//    Input/output, double X[N], the partial or complete candidate vector.
//
//    Input/output, int &INDX, a communication flag.
//    On input,
//      0 to start a search.
//    On output:
//      1, a complete output vector has been determined and returned in X(1:N);
//      2, candidates are needed for position X(K);
//      3, no more possible vectors exist.
//
//    Inout/output, int &K, if INDX=2, the current vector index being considered.
//
//    Input/output, int &NSTACK, the current length of the stack.
//
//    Input/output, int NCAN[N], lists the current number of candidates for
//    positions 1 through K.
//
{
//
//  If this is the first call, request a candidate for position 1.
//
  if ( indx == 0 )
  {
    k = 1;
    nstack = 0;
    indx = 2;
    return;
  }
//
//  Examine the stack.
//
  for ( ; ; )
  {
//
//  If there are candidates for position K, take the first available
//  one off the stack, and increment K.
//
//  This may cause K to reach the desired value of N, in which case
//  we need to signal the user that a complete set of candidates
//  is being returned.
//
    if ( 0 < ncan[k-1] )
    {
      x[k-1] = stack[nstack-1];
      nstack = nstack - 1;

      ncan[k-1] = ncan[k-1] - 1;

      if ( k != n )
      {
        k = k + 1;
        indx = 2;
      }
      else
      {
        indx = 1;
      }

      break;
    }
//
//  If there are no candidates for position K, then decrement K.
//  If K is still positive, repeat the examination of the stack.
//
    else
    {
      k = k - 1;

      if ( k <= 0 )
      {
        indx = 3;
        break;
      }

    }

  }

  return;
}
//****************************************************************************80

double r8vec_frac ( int n, double a[], int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_FRAC searches for the K-th smallest entry in an R8VEC.
//
//  Discussion:
//
//    Hoare's algorithm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input/output, double A[N].
//    On input, A is the array to search.
//    On output, the elements of A have been somewhat rearranged.
//
//    Input, int K, the fractile to be sought.  If K = 1, the minimum
//    entry is sought.  If K = N, the maximum is sought.  Other values
//    of K search for the entry which is K-th in size.  K must be at
//    least 1, and no greater than N.
//
//    Output, double R8VEC_FRAC, the value of the K-th fractile of A.
//
{
  double afrac;
  int i;
  int iryt;
  int j;
  int left;
  double temp;
  double x;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal nonpositive value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( k <= 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal nonpositive value of K = " << k << "\n";
    exit ( 1 );
  }

  if ( n < k )
  {
    cerr << "\n";
    cerr << "R8VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal N < K, K = " << k << "\n";
    exit ( 1 );
  }

  left = 1;
  iryt = n;

  for ( ; ; )
  {
    if ( iryt <= left )
    {
      afrac = a[k-1];
      break;
    }

    x = a[k-1];
    i = left;
    j = iryt;

    for ( ; ; )
    {
      if ( j < i )
      {
        if ( j < k )
        {
          left = i;
        }
        if ( k < i )
        {
          iryt = j;
        }
        break;
      }
//
//  Find I so that X <= A(I).
//
      while ( a[i-1] < x )
      {
        i = i + 1;
      }
//
//  Find J so that A(J) <= X
//
      while ( x < a[j-1] )
      {
        j = j - 1;
      }

      if ( i <= j )
      {
        temp = a[i-1];
        a[i-1] = a[j-1];
        a[j-1] = temp;
        i = i + 1;
        j = j - 1;
      }
    }
  }

  return afrac;
}
//****************************************************************************80

void r8vec_indicator ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDICATOR sets an R8VEC to the indicator vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
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

  for ( i = 0; i < n; i++ ) 
  {
    a[i] = ( double ) ( i + 1 );
  }

  return;
}
//****************************************************************************80

bool r8vec_mirror_next ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
//
//  Discussion:
//
//    In normal use, the user would set every element of A to be positive.
//    The routine will take the input value of A, and output a copy in
//    which the signs of one or more entries have been changed.  Repeatedly
//    calling the routine with the output from the previous call will generate
//    every distinct "variation" of A; that is, all possible sign variations.
//
//    When the output variable DONE is TRUE (or equal to 1), then the
//    output value of A_NEW is the last in the series.
//
//    Note that A may have some zero values.  The routine will essentially
//    ignore such entries; more exactly, it will not stupidly assume that -0
//    is a proper "variation" of 0!
//
//    Also, it is possible to call this routine with the signs of A set
//    in any way you like.  The routine will operate properly, but it
//    will nonethess terminate when it reaches the value of A in which
//    every nonzero entry has negative sign.
//
//
//    More efficient algorithms using the Gray code seem to require internal
//    memory in the routine, which is not one of MATLAB's strong points,
//    or the passing back and forth of a "memory array", or the use of
//    global variables, or unnatural demands on the user.  This form of
//    the routine is about as clean as I can make it.
//
//  Example:
//
//      Input         Output
//    ---------    --------------
//    A            A         R8VEC_MIRROR_NEXT
//    ---------    --------  ----
//     1  2  3     -1  2  3  false
//    -1  2  3      1 -2  3  false
//     1 -2  3     -1 -2  3  false
//    -1 -2  3      1  2 -3  false
//     1  2 -3     -1  2 -3  false
//    -1  2 -3      1 -2 -3  false
//     1 -2 -3     -1 -2 -3  false
//    -1 -2 -3      1  2  3  true
//
//     1  0  3     -1  0  3  false
//    -1  0  3      1  0 -3  false
//     1  0 -3     -1  0 -3  false
//    -1  0 -3      1  0  3  true
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, double A[N], a vector of real numbers.  On output, some
//    signs have been changed.
//
//    Output, bool R8VEC_MIRROR_NEXT, is TRUE if the input vector A was 
//    the last element in the series (every entry was nonpositive); the 
//    output vector is reset so that all entries are nonnegative, but 
//    presumably the ride is over!
//
{
  int i;
  int positive;
//
//  Seek the first strictly positive entry of A.
//
  positive = 0;
  for ( i = 1; i <= n; i++ )
  {
    if ( 0.0 < a[i-1] )
    {
      positive = i;
      break;
    }
  }
//
//  If there is no strictly positive entry of A, there is no successor.
//
  if ( positive == 0 )
  {
    for ( i = 1; i <= n; i++ )
    {
      a[i-1] = -a[i-1];
    }
    return true;
  }
//
//  Otherwise, negate A up to the positive entry.
//
  for ( i = 1; i <= positive; i++ )
  {
    a[i-1] = -a[i-1];
  }
  return false;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
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
  for ( i = 0; i <= n-1; i++ ) 
  {
    cout << "  " << setw(8)  << i 
         << "  " << setw(12) << a[i] << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec_uniform ( int n, double b, double c, int &seed, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM returns a scaled pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2005
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
//    Input, int N, the number of entries in the vector.
//
//    Input, double B, C, the lower and upper limits of the pseudorandom values.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double X[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    x[i] = b + ( c - b ) * ( double ) ( seed ) * 4.656612875E-10;
  }

  return;
}
//****************************************************************************80

void r8vec_uniform_01 ( int n, int &seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2004
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
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return;
}
//****************************************************************************80

unsigned long rand_initialize ( unsigned long seed )

//****************************************************************************80
//
//  Purpose:
//
//    RAND_INITIALIZE initializes the RAND random number generator.
//
//  Discussion:
//
//    If you don't initialize RAND, the random number generator, 
//    it will behave as though it were seeded with value 1.  
//    This routine will either take a user-specified seed, or
//    (if the user passes a 0) make up a "random" one.  In either
//    case, the seed is passed to SRAND (the appropriate routine 
//    to call when setting the seed for RAND).  The seed is also
//    returned to the user as the value of the function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned long SEED, is either 0, which means that the user
//    wants this routine to come up with a seed, or nonzero, in which
//    case the user has supplied the seed.
//
//    Output, unsigned long RAND_INITIALIZE, is the value of the seed
//    passed to SRAND, which is either the user's input value, or if
//    that was zero, the value selected by this routine.
//
{
  if ( seed != 0 )
  {
    cout << "\n";
    cout << "RAND_INITIALIZE\n";
    cout << "  Initialize RAND with user SEED = " << seed << "\n";
  }
  else
  {
    seed = get_seed ( );

    cout << "\n";
    cout << "RAND_INITIALIZE\n";
    cout << "  Initialize RAND with arbitrary SEED = " << seed << "\n";
  }
//
//  Now set the seed.
//
  srand ( seed );

  return seed;
}
//****************************************************************************80

unsigned long random_initialize ( unsigned long seed )

//****************************************************************************80
//
//  Purpose:
//
//    RANDOM_INITIALIZE initializes the RANDOM random number generator.
//
//  Discussion:
//
//    If you don't initialize RANDOM, the random number generator, 
//    it will behave as though it were seeded with value 1.  
//    This routine will either take a user-specified seed, or
//    (if the user passes a 0) make up a "random" one.  In either
//    case, the seed is passed to SRANDOM (the appropriate routine 
//    to call when setting the seed for RANDOM).  The seed is also
//    returned to the user as the value of the function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned long SEED, is either 0, which means that the user
//    wants this routine to come up with a seed, or nonzero, in which
//    case the user has supplied the seed.
//
//    Output, unsigned long RANDOM_INITIALIZE, is the value of the seed
//    passed to SRANDOM, which is either the user's input value, or if
//    that was zero, the value selected by this routine.
//
{
  if ( seed != 0 )
  {
    cout << "\n";
    cout << "RANDOM_INITIALIZE\n";
    cout << "  Initialize RANDOM with user SEED = " << seed << "\n";
  }
  else
  {
    seed = get_seed ( );

    cout << "\n";
    cout << "RANDOM_INITIALIZE\n";
    cout << "  Initialize RANDOM with arbitrary SEED = " << seed << "\n";
  }
//
//  Now set the seed.
//
  srandom ( seed );

  return seed;
}
//****************************************************************************80

void rat_add ( int itop1, int ibot1, int itop2, int ibot2, int &itop, int &ibot, 
  bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_ADD adds two rational values.
//
//  Discussion:
//
//    The routine computes
//
//      ITOP/IBOT = ITOP1/IBOT1 + ITOP2/IBOT2
//
//    while trying to avoid int overflow.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ITOP1, IBOT1, the first value to add.
//
//    Input, int ITOP2, IBOT2, the second value to add.
//
//    Output, int &ITOP, &IBOT, the sum.
//
//    Output, bool &ERROR, is TRUE if an error occurred.
//
{
  int i_max;
  int itemp;
  int ibot3;

  i_max = i4_huge ( );

  error = false;

  if ( itop1 == 0 )
  {
    itop = itop2;
    ibot = ibot2;
    return;
  }
  else if ( itop2 == 0 )
  {
    itop = itop1;
    ibot = ibot1;
    return;
  }
//
//  Divide out  the greatest common factor of the two denominators.
//
  ibot3 = i4_gcd ( ibot1, ibot2 );
  ibot1 = ibot1 / ibot3;
  ibot2 = ibot2 / ibot3;
//
//  The fraction may now be formally written as:
//
//    (itop1*ibot2 + itop2*ibot1) / (ibot1*ibot2*ibot3)
//
//  Check the tops for overflow.
//
  if ( ( double ) ( i_max ) < fabs ( ( double ) ( itop1 ) * ( double ) ( ibot2 ) ) )
  {
    error = true;
    cerr << "\n";
    cerr << "RAT_ADD - Fatal error!\n";
    cerr << "  Overflow of top of rational sum.\n";
    itop = 0;
    ibot = 1;
    exit ( 1 );
  }

  itop1 = itop1 * ibot2;

  if ( ( double ) ( i_max ) < fabs ( ( double ) ( itop2 ) * ( double ) ( ibot1 ) ) )
  {
    error = true;
    cerr << "\n";
    cerr << "RAT_ADD - Fatal error!\n";
    cerr << "  Overflow of top of rational sum.\n";
    itop = 0;
    ibot = 1;
    exit ( 1 );
  }

  itop2 = itop2 * ibot1;

  if ( ( double ) ( i_max ) < fabs ( ( double ) ( itop1 ) + ( double ) ( itop2 ) ) )
  {
    error = true;
    cerr << "\n";
    cerr << "RAT_ADD - Fatal error!\n";
    cerr << "  Overflow of top of rational sum.\n";
    itop = 0;
    ibot = 1;
    exit ( 1 );
  }

  itop = itop1 + itop2;
//
//  Check the bottom for overflow.
//
  if ( ( double ) ( i_max ) < 
    fabs ( ( double ) ( ibot1 ) * ( double ) ( ibot2 ) * ( double ) ( ibot3 ) ) )
  {
    error = true;
    cerr << "\n";
    cerr << "RAT_ADD - Fatal error!\n";
    cerr << "  Overflow of bottom of rational sum.\n";
    itop = 0;
    ibot = 1;
    exit ( 1 );
  }

  ibot = ibot1 * ibot2 * ibot3;
//
//  Put the fraction in lowest terms.
//
  itemp = i4_gcd ( itop, ibot );
  itop = itop / itemp;
  ibot = ibot / itemp;
//
//  The bottom should be positive.
//
  if ( ibot < 0 )
  {
    ibot = -ibot;
    itop = -itop;
  }

  return;
}
//****************************************************************************80

void rat_div ( int itop1, int ibot1, int itop2, int ibot2, int &itop, 
  int &ibot, bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_DIV divides one rational value by another.
//
//  Discussion:
//
//    The routine computes
//
//      ITOP / IBOT = ( ITOP1 / IBOT1 ) / ( ITOP2 / IBOT2 ).
//
//    while avoiding integer overflow.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ITOP1, IBOT1, the numerator.
//
//    Input, int ITOP2, IBOT2, the denominator.
//
//    Output, int &ITOP, &IBOT, the result.
//
//    Output, bool &ERROR, is TRUE if an error occurred.
//
{
  int i_max;
  int itemp;

  error = false;

  i_max = i4_huge ( );

  if ( ibot1 == 0 || itop2 == 0 || ibot2 == 0 )
  {
    error = true;
    return;
  }

  if ( itop1 == 0 )
  {
    itop = 0;
    ibot = 1;
    return;
  }
//
//  Get rid of all common factors in top and bottom.
//
  itemp = i4_gcd ( itop1, ibot1 );
  itop1 = itop1 / itemp;
  ibot1 = ibot1 / itemp;
  itemp = i4_gcd ( itop1, itop2 );
  itop1 = itop1 / itemp;
  itop2 = itop2 / itemp;
  itemp = i4_gcd ( ibot2, ibot1 );
  ibot2 = ibot2 / itemp;
  ibot1 = ibot1 / itemp;
  itemp = i4_gcd ( ibot2, itop2 );
  ibot2 = ibot2 / itemp;
  itop2 = itop2 / itemp;
//
//  The fraction (ITOP1*IBOT2)/(IBOT1*ITOP2) is in lowest terms.
//
//  Check the top for overflow.
//
  if ( ( double ) i_max < fabs ( ( double ) itop1 * ( double ) ibot2 ) )
  {
    error = true;
    cerr << "\n";
    cerr << "RAT_DIV - Fatal error!\n";
    cerr << "  Overflow of top of rational fraction.\n";
    itop = 0;
    ibot = 0;
    exit ( 1 );
  }

  itop = itop1 * ibot2;
//
//  Check the bottom IBOT1*ITOP2 for overflow.
//
  if ( ( double ) i_max < fabs ( ( double ) ibot1 * ( double ) itop2 ) )
  {
    error = true;
    cerr << "\n";
    cerr << "RAT_DIV - Fatal error!\n";
    cerr << "  Overflow of bottom of rational fraction.\n";
    itop = 0;
    ibot = 1;
    exit ( 1 );
  }
  ibot = ibot1 * itop2;
//
//  The bottom should be positive.
//
  if ( ibot < 0 )
  {
    ibot = - ibot;
    itop = - itop;
  }
//
//  The fraction is ITOP/IBOT with no loss of accuracy.
//
  return;
}
//****************************************************************************80

void rat_farey ( int n, int max_frac, int &num_frac, int a[], int b[] )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_FAREY computes the N-th row of the Farey fraction table.
//
//  Example:
//
//    N = 5
//
//    NUM_FRAC = 11
//    A =  0  1  1  1  2  1  3  2  3  4  1
//    B =  1  5  4  3  5  2  5  3  4  5  1
//
//  Discussion:
//
//    In this form of the Farey fraction table, fractions in row N lie between
//    0 and 1, are in lowest terms, and have a denominator that is no greater
//    than N.  Row N is computed directly, and does not require the computation
//    of previous rows.
//
//    The data satisfy the relationship:
//
//      A(K+1) * B(K) - A(K) * B(K+1) = 1
//
//    The number of items in the N-th row is roughly N**2 / PI**2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 1, Fundamental Algorithms,
//    Addison Wesley, 1968, page 157.
//
//  Parameters:
//
//    Input, int N, the desired row number.  N must be positive.
//
//    Input, int MAX_FRAC, the maximum number of fractions to compute.
//
//    Output, int &NUM_FRAC, the number of fractions computed.
//
//    Output, int A[MAX_FRAC], B[MAX_FRAC], contains the NUM_FRAC
//    numerators and denominators of the N-th row of the Farey fraction table.
//
{
  int c;
  int k;

  if ( n <= 0 )
  {
    num_frac = 0;
    return;
  }

  k = 0;

  if ( max_frac <= 0 )
  {
    num_frac = k;
    return;
  }

  a[k] = 0;
  b[k] = 1;
  k = 1;

  if ( max_frac <= 1 )
  {
    num_frac = k;
    return;
  }

  a[k] = 1;
  b[k] = n;
  k = 2;

  while ( k < max_frac )
  {
    if ( a[k-1] == 1 && b[k-1] == 1 )
    {
      break;
    }

    c = ( b[k-2] + n ) / b[k-1];
    a[k] = c * a[k-1] - a[k-2];
    b[k] = c * b[k-1] - b[k-2];
    k = k + 1;
  }

  num_frac = k;

  return;
}
//****************************************************************************80

void rat_farey2 ( int n, int a[], int b[] )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_FAREY2 computes the next row of the Farey fraction table.
//
//  Example:
//
//    Input:
//
//      N = 3
//      A =  0  1  1  2  1
//      B =  1  3  2  3  1
//
//    Output:
//
//      A =  0  1  1  2  1  3  2  3  1
//      B =  1  4  3  5  2  5  3  4  1
//
//  Discussion:
//
//    In this form of the Farey fraction table, fractions in row N lie between
//    0 and 1, and are in lowest terms.  For every adjacent pair of input
//    fractions, A1/B1 and A2/B2, the mediant (A1+A2)/(B1+B2) is computed
//    and inserted between them.
//
//    The number of items in the N-th row is 1+2**(N-1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the input row number.  N must be nonnegative.
//    If N is zero, then the input is ignored, and the entries of
//    row 1 are computed directly.
//
//    Input/output, int A[1+2**N], B[1+2**N].
//    On input, entries 1 through 1+2**(N-1) contain the entries of row N.
//    On output, entries 1 through 1+2**N contain the entries of row N+1.
//
{
  int i;
  int ihi;

  if ( n == 0 )
  {
    a[0] = 0;
    b[0] = 1;
    a[1] = 1;
    b[1] = 1;
    return;
  }
//
//  Shift the current data.
//
  ihi = ( int ) pow ( ( double ) 2, n-1 );
  for ( i = ihi; 0 <= i; i-- )
  {
    a[2*i] = a[i];
    b[2*i] = b[i];
  }
//
//  Compute the mediants.
//
  ihi = ( int ) pow ( ( double ) 2, n ) - 1;

  for ( i = 1; i <= ihi; i = i + 2 )
  {
    a[i] = a[i-1] + a[i+1];
    b[i] = b[i-1] + b[i+1];
  }

  return;
}
//****************************************************************************80

void rat_mul ( int itop1, int ibot1, int itop2, int ibot2, int &itop, 
  int &ibot, bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_MUL multiplies two fractions.
//
//  Discussion:
//
//    The routine computes
//
//      ITOP / IBOT = ( ITOP1 / IBOT1 ) * ( ITOP2 / IBOT2 ).
//
//    while avoiding int overflow.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ITOP1, IBOT1, the first factor.
//
//    Input, int ITOP2, IBOT2, the second factor.
//
//    Output, int &ITOP, &IBOT, the product.
//
//    Output, bool &ERROR, is TRUE if an error occurred.
//
{
  int i_max;
  int itemp;

  error = false;

  i_max = i4_huge ( );

  if ( itop1 == 0 || itop2 == 0 )
  {
    itop = 0;
    ibot = 1;
    return;
  }
//
//  Get rid of all common factors in top and bottom.
//
  itemp = i4_gcd ( itop1, ibot1 );
  itop1 = itop1 / itemp;
  ibot1 = ibot1 / itemp;
  itemp = i4_gcd ( itop1, ibot2 );
  itop1 = itop1 / itemp;
  ibot2 = ibot2 / itemp;
  itemp = i4_gcd ( itop2, ibot1 );
  itop2 = itop2 / itemp;
  ibot1 = ibot1 / itemp;
  itemp = i4_gcd ( itop2, ibot2 );
  itop2 = itop2 / itemp;
  ibot2 = ibot2 / itemp;
//
//  The fraction (ITOP1 * ITOP2) / (IBOT1*IBOT2) is in lowest terms.
//
//  Check the top ITOP1 * ITOP2 for overflow.
//
  if ( ( double ) ( i_max ) < fabs ( ( double ) itop1 * ( double ) itop2 ) )
  {
    error = true;
    cerr << "\n";
    cerr << "RAT_MUL - Fatal error!\n";
    cerr << "  Overflow of top of rational product.\n";
    itop = 0;
    exit ( 1 );
  }

  itop = itop1 * itop2;
//
//  Check the bottom IBOT1*IBOT2 for overflow.
//
  if ( ( double ) ( i_max ) < fabs ( ( double ) ibot1 * ( double ) ibot2 ) )
  {
    error = true;
    cerr << "\n";
    cerr << "RAT_MUL - Fatal error!\n";
    cerr << "  Overflow of bottom of rational product.\n";
    ibot = 1;
    exit ( 1 );
  }

  ibot = ibot1 * ibot2;
//
//  The bottom should be positive.
//
  if ( ibot < 0 )
  {
    ibot = -ibot;
    itop = -itop;
  }

  return;
}
//****************************************************************************80

void rat_normalize ( int &a, int &b )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_NORMALIZE normalizes a rational number.
//
//  Discussion:
//
//    If A = B = 0, return.
//
//    If A = 0 (and B nonzero) set B => 1 and return.
//
//    If A nonzero, and B = 0, then A => 1 and return.
//
//    If A = B, then set A => 1, B => 1 and return.
//
//    If B < 0, then A => -A, B => -B.
//
//    If 1 < C = GCD(|A|,|B|), A => A/C, B => B/C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &A, &B, the rational number.
//
{
  int c;
//
//  Cases where one or both is 0.
//
  if ( a == 0 && b == 0 )
  {
    return;
  }
  else if ( a == 0 && b != 0 )
  {
    b = 1;
    return;
  }
  else if ( a != 0 && b == 0 )
  {
    a = 1;
    return;
  }

  if ( a == b )
  {
    a = 1;
    b = 1;
    return;
  }

  if ( b < 0 )
  {
    a = -a;
    b = -b;
  }

  c = i4_gcd ( abs ( a ), abs ( b ) );

  if ( 1 < c )
  {
    a = a / c;
    b = b / c;
  }

  return;
}
//****************************************************************************80

void rat_sum_formula ( int n, int a[], int b[] )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_SUM_FORMULA computes the formulas for sums of powers of integers.
//
//  Example:
//
//    N = 6
//
//        1    2    3    4    5    6    7
//    -----------------------------------
//    0 | 1    0    0    0    0    0    0
//      |
//    1 | 1    1    0    0    0    0    0
//      | 2    2
//      |
//    2 | 1    1    1    0    0    0    0
//      | 3    2    6
//      |
//    3 | 1    1    1    0    0    0    0
//      | 4    2    4
//      | 
//    4 | 1    1    1    0   -1    0    0
//      | 5    2    3        30
//      |
//    5 | 1    1    5    0   -1    0    0
//      | 6    2   12        12
//      |
//    6 | 1    1    1    0   -1    0    1
//      | 7    2    2         6        42
//
//    The interpretation of row 2, for instance, is:
//
//      sum ( 1 <= I <= N ) I**2 = 1/3 N**3 + 1/2 N**2 + 1/6 N
//
//    This suggests that a more sensible way to display the table
//    is to reverse the order of the entries in the row, so that
//    the entry in column J is the coeficient of N**J, which is
//    not the case now.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Owens,
//    Sums of Powers of Integers,
//    Mathematics Magazine,
//    Volume 65, Number 1, February 1992, pages 38-40.
//
//  Parameters:
//
//    Input, int N, the number of rows of coefficients to compute.
//
//    Output, int A[(N+1)*(N+1)], B[(N+1)*(N+1)], the numerator and denominator
//    of the coefficients.
//
{
  int asum;
  int aval;
  int bsum;
  int bval;
  bool error;
  int i;
  int j;

  a[0+0*(n+1)] = 1;
  for ( j = 1; j < n+1; j++ )
  {
    a[0+j*(n+1)] = 0;
  }

  b[0+0*(n+1)] = 1;
  for ( j = 1; j < n+1; j++ )
  {
    b[0+j*(n+1)] = 1;
  }

  for ( i = 1; i <= n; i++ )
  {
    asum = 0;
    bsum = 0;
//
//  Subdiagonal entries are multiples of entries above them.
//
    for ( j = 0; j < i; j++ )
    {
      aval = a[i-1+j*(n+1)];
      bval = b[i-1+j*(n+1)];

      rat_mul ( aval, bval, i, i+1-j, aval, bval, error );

      a[i+j*(n+1)] = aval;
      b[i+j*(n+1)] = bval;

      rat_add ( asum, bsum, aval, bval, asum, bsum, error );
    }
//
//  Diagonal entry is 1 - sum of previous entries in row.
//
    asum = -asum;

    rat_add ( 1, 1, asum, bsum, aval, bval, error );

    a[i+i*(n+1)] = aval;
    b[i+i*(n+1)] = bval;
//
//  Superdiagonal entries are zero.
//
    for ( j = i+1; j < n+1; j++ )
    { 
      a[i+j*(n+1)] = 0;
    }

    for ( j = i+1; j <= n; j++ )
    { 
      b[i+j*(n+1)] = 1;
    }
  }

  return;
}
//****************************************************************************80

void rat_to_cfrac ( int ip, int iq, int m, int &n, int a[], bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_TO_CFRAC converts a rational value to a continued fraction.
//
//  Discussion:
//
//    The routine is given a rational number represented by IP/IQ, and
//    computes the monic or "simple" continued fraction representation
//    with int coefficients of the number:
//
//      A(1) + 1/ (A(2) + 1/ (A(3) + ... + 1/A(N) ...))
//
//    The user must dimension A to a value M which is "large enough".
//    The actual number of terms needed in the continued fraction
//    representation cannot be known beforehand.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
//    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, 
//    Christoph Witzgall.
//    C++ version by John Burkardt
//
//  Reference:
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
//    John Rice, Henry Thatcher, Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968.
//
//  Parameters:
//
//    Input, int IP, IQ, the numerator and denominator of the
//    rational value whose continued fraction representation is
//    desired.
//
//    Input, int M, the dimension of A.  If M is not great
//    enough, the algorithm may run out of space.
//
//    Output, int &N, the actual number of entries used in A.
//
//    Output, int A[M], contains the continued fraction
//    representation of the number.
//
//    Output, bool &ERROR, is TRUE if an error occurred.
//
{
  error = false;

  n = 0;

  for ( ; ; )
  {
    if ( m <= n )
    {
      error = true;
      return;
    }

    a[n] = ip / iq;
    n = n + 1;

    ip = ip % iq;

    if ( ip == 0 )
    {
      return;
    }

    if ( m <= n )
    {
      error = true;
      return;
    }

    a[n] = iq / ip;
    n = n + 1;

    iq = iq % ip;

    if ( iq == 0 )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void rat_to_dec ( int rat_top, int rat_bot, int &mantissa, int &exponent )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_TO_DEC converts a rational to a decimal representation.
//
//  Discussion:
//
//    A rational value is represented by RAT_TOP / RAT_BOT.
//
//    A decimal value is represented as MANTISSA * 10**EXPONENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int RAT_TOP, RAT_BOT, the rational value.
//
//    Output, int &MANTISSA, &EXPONENT, the decimal number.
//
{
  int gcd;
  double r;
  double r_max;
  int s;
//
//  Handle an input value of 0.
//
  if ( rat_top == 0 ) 
  {
    mantissa = 0;
    exponent = 0;
    return;
  }
//
//  Take out any common factor.
//
  gcd = i4_gcd ( rat_top, rat_bot );
  rat_top = rat_top / gcd;
  rat_bot = rat_bot / gcd;
//
//  Force the bottom of the fraction to be positive.
//
  if ( rat_bot < 0 )
  {
    rat_top = -rat_top;
    rat_bot = -rat_bot;
  }
//
//  If the fraction is a whole number, we can finish now.
//
  if ( rat_bot == 1 )
  {
    mantissa = rat_top;
    exponent = 0;
    return;
  }
//
//  Whole factors of 10 in the bottom or top go into the decimal exponent.
//
  exponent = 0;

  while ( ( rat_bot % 10 ) == 0 )
  {
    exponent = exponent - 1;
    rat_bot = rat_bot / 10;
  }

  while ( ( rat_top % 10 ) == 0 )
  {
    exponent = exponent + 1;
    rat_top = rat_top / 10;
  }
//
//  Convert to a real value.
//
  r = ( double ) rat_top / ( double ) rat_bot;

  if ( r < 0 )
  {
    s = -1;
    r = -r;
  }
  else
  {
    s = 1;
  }

  r_max = ( ( double ) i4_huge ( ) ) / 10.0;

  while ( r != ( double ) ( ( int ) r ) && r < r_max )
  {
    r = r * 10.0;
    exponent = exponent - 1;
  }

  mantissa = s * ( int ) r;

  return;
}
//****************************************************************************80

double rat_to_r8 ( int top, int bot )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_TO_R8 converts rational values to R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TOP, BOT, the rational quantity
//    (TOP/BOT) that is to be converted.
//
//    Output, double RAT_TO_R8, the real value of the fraction.
//
{
  return ( ( double ) top / ( double ) bot );
}
//****************************************************************************80

int rat_width ( int a, int b )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_WIDTH returns the "width" of a rational number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, B, the rational number.
//
//    Output, int RAT_WIDTH, the "width" of the rational number.
//
{
  int abs_a;
  int abs_b;
  int ten_pow;
  int value;

  value = 1;
  ten_pow = 10;

  if ( a == 0 )
  {
    return value;
  }
  
  abs_a = abs ( a );

  while ( ten_pow <= abs_a )
  {
    value = value + 1;
    ten_pow = ten_pow * 10;
  }
//
//  If the fraction is negative, a minus sign will be prepended to the
//  numerator.
//
  if ( a * b < 0 )
  {
    value = value + 1;
    ten_pow = ten_pow * 10;
  }

  abs_b = abs ( b );

  while ( ten_pow <= abs_b )
  {
    value = value + 1;
    ten_pow = ten_pow * 10;
  }

  return value;
}
//****************************************************************************80

void ratmat_det ( int n, int iatop[], int iabot[], int &idtop, int &idbot, 
  bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    RATMAT_DET finds the determinant of an N by N matrix of rational entries.
//
//  Discussion:
//
//    The brute force method is used.
//
//    This routine should only be used for small matrices, since this
//    calculation requires the summation of N! products of N numbers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of A.
//
//    Input, int IATOP[N*N], IABOT[N*N], the numerators
//    and denominators of the entries of the matrix.
//
//    Output, int &IDTOP, &IDBOT, the determinant of the matrix,
//    expressed as IDTOP/IDBOT.
//
//    Output, bool &ERROR, is TRUE if an error occurred.
//
{
  bool even;
  int i;
  int *iarray;
  int ibot;
  int ibot2;
  int itop;
  int itop2;
  bool more;

  error = false;

  more = false;
  idtop = 0;
  idbot = 1;

  iarray = new int[n];

  for ( ; ; )
  {
    perm_next ( n, iarray, more, even );

    if ( even )
    {
      itop = 1;
    }
    else
    {
      itop = -1;
    }

    ibot = 1;

    for ( i = 1; i <= n; i++ )
    {
      itop2 = iatop[i-1+(iarray[i-1]-1)*n];
      ibot2 = iabot[i-1+(iarray[i-1]-1)*n];

      rat_mul ( itop, ibot, itop2, ibot2, itop, ibot, error );

      if ( error )
      {
        cerr << "\n";
        cerr << "RATMAT_DET - Fatal error!\n";
        cerr << "  An overflow occurred.\n";
        cerr << "  The determinant calculation cannot be done\n";
        cerr << "  for this matrix.\n";
        idtop = 0;
        idbot = 1;
        delete [] iarray;
        exit ( 1 );
      }
    }

    rat_add ( itop, ibot, idtop, idbot, idtop, idbot, error );

    if ( error )
    {
      cerr << "\n";
      cerr << "RATMAT_DET - Fatal error!\n";
      cerr << "  An overflow occurred.\n";
      cerr << "  The determinant calculation cannot be done\n";
      cerr << "  for this matrix.\n";
      idtop = 0;
      idbot = 1;
      delete [] iarray;
      exit ( 1 );
    }

    if ( !more )
    {
      break;
    }
  }

  delete [] iarray;
//
//  The bottom should be positive.
//
  if ( idbot < 0 )
  {
    idbot = -idbot;
    idtop = -idtop;
  }

  return;
}
//****************************************************************************80

void ratmat_print ( int m, int n, int a[], int b[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    RATMAT_PRINT prints out a rational matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, int A[M*N], B[M*N], the rational matrix.
//
//    Input, string TITLE, a title.
//
{
  bool all_ones;
  int i;
  int j;
  int jmax;
  int jmin;
  int kmax;
  int ncolum = 80;
  int npline;
//
//  Figure out how many rationals we can get in NCOLUM columns.
//
//  Can I pass the VALUE of these entries this way?
//
  kmax = 0;

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      kmax = i4_max ( kmax, rat_width ( a[i+j*m], b[i+j*m] ) );
    }
  }
  kmax = kmax + 2;
  npline = ncolum / kmax;
//
//  Now do the printing.
//
  for ( jmin = 0; jmin < n; jmin = jmin + npline )
  {
    jmax = i4_min ( jmin+npline-1, n-1 );

    cout << "\n";

    if ( jmin == 0 )
    {
      cout << title << "\n";
      cout << "\n";
    }

    if ( 0 < jmin || jmax < n-1 )
    {
      cout << "Columns " << jmin << " to " << jmax << "\n";
      cout << "\n";
    }

    for ( i = 0; i < m; i++ )
    {
      for ( j = jmin; j <= jmax; j++ )
      {
        cout << setw(kmax) << a[i+j*m] << "  ";
      }
      cout << "\n";
//
//  Delete each denominator that is 1.  If all are 1, don't
//  even print out the line.
//
      all_ones = true;

      for ( j = jmin; j <= jmax; j++ )
      {
        if ( b[i+j*m] != 1 )
        {
          all_ones = false;
        }
      }

      if ( !all_ones )
      {
        for ( j = jmin; j <= jmax; j++ )
        {
          cout << setw(kmax) << b[i+j*m] << "  ";
        }
        cout << "\n";
      }

      if ( jmax == n - 1 && i == m - 1 )
      {
      }
      else
      {
        cout << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void regro_next ( bool &done, int n, int v[], int vmax[] )

//****************************************************************************80
//
//  Purpose:
//
//    REGRO_NEXT computes restricted growth functions one at a time.
//
//  Discussion:
//
//    A restricted growth function on N is a vector (V(1), ..., V(N) )
//    of values V(I) between 1 and N, satisfying the requirements:
//      V(1) = 1;
//      V(I) <= 1 + max ( V(1), V(2), ..., V(I-1) ).
//
//    The number of restricted growth functions on N is equal to
//    the Bell number B(N).
//
//    There is a bijection between restricted growth functions on N
//    and set partitions of N.
//
//  Example:
//
//    The 15 restricted growth functions for N = 4 are:
//
//    (1111), (1112), (1121), (1122), (1123),
//    (1211), (1212), (1213), (1221), (1222),
//    (1223), (1231), (1232), (1233), (1234).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input/output, bool &DONE.
//    On first call, set DONE to TRUE, and then do not alter it.
//    On output, DONE will be FALSE if the routine has computed another
//    restricted growth function, or TRUE if all the restricted
//    growth functions have been returned.
//
//    Input, int N, the number of components in the restricted growth
//    function.
//
//    Input/output, int V[N].  The user need not set this quantity
//    before the initial call, and should not alter it between successive
//    calls.  On each return from the routine, with DONE = .FALSE.,
//    V will contain the componentwise values of the next restricted
//    growth function.
//
//    Input/output, int VMAX[N].  The user need not set this quantity
//    before the initial call, and should not alter it between calls.
//    VMAX(I) records the largest value that component V(I) could take,
//    given the values of components 1 through I-1.
//
{
  int i;
  int j;
//
//  First call:
//
  if ( done )
  {
    for ( i = 0; i < n; i++ )
    {
      v[i] = 1;
    }

    vmax[0] = 1;
    for ( i = 1; i < n; i++ )
    {
      vmax[i] = 2;
    }

    done = false;
  }
//
//  Later calls.
//
  else
  {
    j = n;

    for ( ; ; )
    {
      if ( j == 1 )
      {
        done = true;
        return;
      }
      if ( v[j-1] != vmax[j-1] )
      {
        break;
      }

      j = j - 1;

    }

    v[j-1] = v[j-1] + 1;

    for ( i = j + 1; i <= n; i++ )
    {
      v[i-1] = 1;

      if ( v[j-1] == vmax[j-1] )
      {
        vmax[i-1] = vmax[j-1] + 1;
      }
      else
      {
        vmax[i-1] = vmax[j-1];
      }
    }
  }

  return;
}
//****************************************************************************80

void rfrac_to_cfrac ( int m, double p[], double q[], double t[], bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    RFRAC_TO_CFRAC converts a rational polynomial fraction to a continued fraction.
//
//  Discussion:
//
//    That is, it accepts
//
//      P(1) + P(2) * X + ... + P(M) * X**(M-1)
//      -------------------------------------------------------
//      Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
//
//    and returns the equivalent continued fraction:
//
//      1 / (T(1) + X/(T(2) + X/(...T(2*M-1)+X/(T(2*M) ... )))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 June 2003
//
//  Author:
//
//    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
//    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
//    John Rice, Henry Thatcher, Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968.
//
//  Parameters:
//
//    Input, int M, defines the number of P coefficients,
//    and is one less than the number of Q coefficients, and one
//    half the number of T coefficients.
//
//    Input, double P[M], Q[M+1], the coefficients defining the rational
//    polynomial fraction.
//
//    Output, double T[2*M], the coefficients defining the continued fraction.
//
//    Output, bool &ERROR, is TRUE if an error occurred.
//
{
  double *a;
  int i;
  int k;
  double ta;

  a = new double[(m+1)*(2*m+1)];

  error = false;

  for ( i = 0; i <= m; i++ )
  {
    a[i+0*(m+1)] = q[i];
  }

  for ( i = 0; i < m; i++ )
  {
    a[i+1*(m+1)] = p[i];
  }

  t[0] = a[0+0*(m+1)] / a[0+1*(m+1)];
  ta = a[m+0*(m+1)];

  for ( i = 1; i <= m; i++ )
  {
    a[m-i+(2*i)*(m+1)] = ta;
  }

  for ( k = 1; k <= 2*m-2; k++ )
  {
    for ( i = 1; i <= (2*m-k)/2; i++ )
    {
      a[i-1+(k+1)*(m+1)] = a[i+(k-1)*(m+1)] - t[k-1] * a[i+k*(m+1)];
    }

    if ( a[0+(k+1)*(m+1)] == 0.0 )
    {
      error = true;
      cerr << "\n";
      cerr << "RFRAC_TO_CFRAC - Fatal error!\n";
      cerr << "  A[0,K+1] is zero for K = " << k << "\n";
      exit ( 1 );
    }

    t[k] = a[0+k*(m+1)] / a[0+(k+1)*(m+1)];
  }

  t[2*m-1] = a[0+(2*m-1)*(m+1)] / a[0+(2*m)*(m+1)];

  delete [] a;

  return;
}
//****************************************************************************80

void rfrac_to_jfrac ( int m, double p[], double q[], double r[], double s[] )

//****************************************************************************80
//
//  Purpose:
//
//    RFRAC_TO_JFRAC converts a rational polynomial fraction to a J fraction.
//
//  Discussion:
//
//    The routine accepts
//
//    P(1) + P(2) * X + ... + P(M) * X**(M-1)
//    -------------------------------------------------------
//    Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
//
//    and returns the equivalent J-fraction:
//
//    R(1) / ( X + S(1) + 
//    R(2) / ( X + S(2) + 
//    R(3) / ...        +
//    R(M) / ( X + S(M) )... ))
//
//    Thanks to Henry Amuasi for noticing and correcting an error in a
//    previous formulation of this routine, 02 October 2010.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
//    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
//    John Rice, Henry Thatcher, Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968.
//
//  Parameters:
//
//    Input, int M, defines the number of P, R, and S coefficients,
//    and is one less than the number of Q coefficients.
//    1 <= M.
//
//    Input, double P[M], Q[M+1], the coefficients defining the rational
//    polynomial fraction.
//
//    Output, double R[M], S[M], the coefficients defining the
//    J-fraction.
//
{
  double *a;
  int i;
  int k;

  if ( m < 1 )
  {
    cerr << "\n";
    cerr << "RFRAC_TO_JFRAC - Fatal error!\n";
    cerr << "  Input M < 1.\n";
    exit ( 1 );
  }

  a = new double[(m+1)*(m+1)];

  for ( i = 0; i <= m; i++ )
  {
    a[i+0*(m+1)] = q[i];
  }

  for ( i = 0; i < m; i++ )
  {
    a[i+1*(m+1)] = p[i];
  }

  if ( 1 < m )
  {
    r[0] = a[m-1+1*(m+1)] / a[m+0*(m+1)];
    s[0] = ( r[0] * a[m-1+0*(m+1)] - a[m-2+1*(m+1)] ) / a[m-1+1*(m+1)];

    for ( k = 0; k < m - 2; k++ )
    {
      a[0+(k+2)*(m+1)] = r[k] * a[0+k*(m+1)] - s[k] * a[0+(k+1)*(m+1)];

      for ( i = 1; i < m - k; i++ )
      {
        a[i+(k+2)*(m+1)] = r[k] * a[i+k*(m+1)] 
          - a[i-1+(k+1)*(m+1)] - s[k] * a[i+(k+1)*(m+1)];
      }

      if ( a[m-2-k+(k+2)*(m+1)] == 0.0 )
      {
        cerr << "\n";
        cerr << "RFRAC_TO_JFRAC - Fatal error!\n";
        cerr << "  A[M-2-K,K+2] = 0 for K = " << k << "\n";
        exit ( 1 );
      }

      r[k+1] = a[m-k-2+(k+2)*(m+1)] / a[m-2-k+1+(k+1)*(m+1)];
      s[k+1] = ( r[k+1] * a[m-2-k+(k+1)*(m+1)] 
        - a[m-2-k-1+(k+2)*(m+1)] ) / a[m-2-k+(k+2)*(m+1)];
    }
    a[0+m*(m+1)] = r[m-2] * a[0+(m-2)*(m+1)] - s[m-2] * a[0+(m-1)*(m+1)];
  }

  r[m-1] = a[0+m*(m+1)] / a[1+(m-1)*(m+1)];
  s[m-1] = a[0+(m-1)*(m+1)] / a[1+(m-1)*(m+1)];

  delete [] a;

  return;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

void schroeder ( int n, int s[] )

//****************************************************************************80
//
//  Purpose:
//
//    SCHROEDER generates the Schroeder numbers.
//
//  Discussion:
//
//    The Schroeder number S(N) counts the number of ways to insert
//    parentheses into an expression of N items, with two or more items within
//    a parenthesis.
//
//    Note that the Catalan number C(N) counts the number of ways
//    to legally arrange a set of N left and N right parentheses.
//
//    The formula is:
//
//      S(N) = ( P(N)(3.0) - 3 P(N-1)(3.0) ) / ( 4 * ( N - 1 ) )
//
//    where P(N)(X) is the N-th Legendre polynomial.
//
//  Example:
//
//    N = 4
//
//    1234
//    12(34)
//    1(234)
//    1(2(34))
//    1(23)4
//    1((23)4)
//    (123)4
//    (12)34
//    (12)(34)
//    (1(23))4
//    ((12)3)4
//
//  First Values:
//
//           1
//           1
//           3
//          11
//          45
//         197
//         903
//        4279
//       20793
//      103049
//      518859
//     2646723
//    13648869
//    71039373
//
//  Recursion:
//
//    S(1) = 1
//    S(2) = 1
//    S(N) = ( ( 6 * N - 9 ) * S(N-1) - ( N - 3 ) * S(N-2) ) / N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    R P Stanley,
//    Hipparchus, Plutarch, Schroeder, and Hough,
//    American Mathematical Monthly,
//    Volume 104, Number 4, 1997, pages 344-350.
//
//    Laurent Habsieger, Maxim Kazarian, Sergei Lando,
//    On the Second Number of Plutarch,
//    American Mathematical Monthly, May 1998, page 446.
//
//  Parameters:
//
//    Input, int N, the number of Schroeder numbers desired.
//
//    Output, int S[N], the Schroeder numbers.
//
{
  int i;

  if ( n <= 0 )
  {
    return;
  }

  s[0] = 1;

  if ( n <= 1 )
  {
    return;
  }

  s[1] = 1;

  if ( n <= 2 )
  {
    return;
  }

  for ( i = 3; i <= n; i++ )
  {
    s[i-1] = ( ( 6 * i - 9 ) * s[i-2] - ( i - 3 ) * s[i-3] ) / i;
  }

  return;
}
//****************************************************************************80

void sort_heap_external ( int n, int &indx, int &i, int &j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int &INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int &I, &J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( indx == 0 )
  {
    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( indx < 0 )
  {
    if ( indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      indx = 2;
      i = i_save;
      j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        indx = 1;
      }
      i = i_save;
      j = j_save;
      return;
    }

    k = k - 1;
    k1 = k;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( indx == 1 ) 
  {
    k1 = k;
  }

  for ( ;; )
  {
    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      indx = -2;
      i = i_save;
      j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    indx = 0;
    i = i_save;
    j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    indx = 1;
    i = i_save;
    j = j_save;
  }

  return;
}
//****************************************************************************80

void subcomp_next ( int n, int k, int a[], bool &more, int &h, int &t )

//****************************************************************************80
//
//  Purpose:
//
//    SUBCOMP_NEXT computes the next subcomposition of N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to a value of N.
//
//    A subcomposition of the integer N into K parts is a composition
//    of M into K parts, where 0 <= M <= N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2008
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer whose subcompositions are desired.
//
//    Input, int K, the number of parts in the subcomposition.
//
//    Input/output, int A[K], the parts of the subcomposition.
//
//    Input/output, bool &MORE, set by the user to start the computation,
//    and by the routine to terminate it.
//
//    Input/output, int &H, &T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;
  static bool more2 = false;
  static int n2 = 0;
//
//  The first computation.
//
  if ( !more )
  {
    n2 = 0;

    for ( i = 0; i < k; i++ )
    {
      a[i] = 0;
    }
    more2 = false;
    h = 0;
    t = 0;

    more = true;
  }
//
//  Do the next element at the current value of N.
//
  else if ( more2 )
  {
    comp_next ( n2, k, a, more2, h, t );
  }
  else
  {
    more2 = false;
    n2 = n2 + 1;

    comp_next ( n2, k, a, more2, h, t );
  }
//
//  Termination occurs if MORE2 = FALSE and N2 = N.
//
  if ( !more2 && n2 == n )
  {
    more = false;
  }

  return;
}
//****************************************************************************80

void subcompnz_next ( int n, int k, int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    SUBCOMPNZ_NEXT computes the next subcomposition of N into K nonzero parts.
//
//  Discussion:
//
//    A composition of the integer N into K nonzero parts is an ordered sequence
//    of K positive integers which sum to a value of N.
//
//    A subcomposition of the integer N into K nonzero parts is a composition
//    of M into K nonzero parts, where 0 < M <= N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer whose subcompositions are desired.
//
//    Input, int K, the number of parts in the subcomposition.
//    K must be no greater than N.
//
//    Input/output, int A[K], the parts of the subcomposition.
//
//    Input/output, bool &MORE, set by the user to start the computation,
//    and by the routine to terminate it.
//
{
  int i;
  static bool more2 = false;
  static int n2 = 0;

  if ( n < k )
  {
    for ( i = 0; i < k; i++ )
    {
      a[i] = -1;
    }
    return;
  }
//
//  The first computation.
//
  if ( !more )
  {
    more = true;
    for ( i = 0; i < k; i++ )
    {
      a[i] = 1;
    }
    n2 = k;
    more2 = false;
  }
//
//  Do the next element at the current value of N.
//
  else if ( more2 )
  {
    compnz_next ( n2, k, a, more2 );
  }
  else
  {
    more2 = false;
    n2 = n2 + 1;

    compnz_next ( n2, k, a, more2 );
  }
//
//  Termination occurs if MORE2 = FALSE and N2 = N.
//
  if ( !more2 && n2 == n )
  {
    more = false;
  }
  return;
}
//****************************************************************************80

void subcompnz2_next ( int n_lo, int n_hi, int k, int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    SUBCOMPNZ2_NEXT computes the next subcomposition of N into K nonzero parts.
//
//  Discussion:
//
//    A composition of the integer N into K nonzero parts is an ordered sequence
//    of K positive integers which sum to a value of N.
//
//    A subcomposition of the integer N into K nonzero parts is a composition
//    of M into K nonzero parts, where 0 < M <= N.
//
//    This routine computes all compositions of K into nonzero parts which sum
//    to values between N_LO and N_HI.
//
//    The routine SUBCOMPNZ_NEXT can be regarded as a special case where N_LO = K.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N_LO, N_HI, the range of values of N for which compositions are desired.
//    N_LO must be no greater than N_HI.
//
//    Input, int K, the number of parts in the subcomposition.
//    K must be no greater than N_HI.
//
//    Input/output, int A[K], the parts of the subcomposition.
//
//    Input/output, bool &MORE, set by the user to start the computation,
//    and by the routine to terminate it.
//
{
  int i;
  static bool more2 = false;
  static int n2 = 0;

  if ( n_hi < k )
  {
    for ( i = 0; i < k; i++ )
    {
      a[i] = -1;
    }
    return;
  }

  if ( n_hi < n_lo )
  {
    for ( i = 0; i < k; i++ )
    {
      a[i] = -1;
    }
    return;
  }
//
//  The first computation.
//
  if ( !more )
  {
    more = true;

    n2 = i4_max ( k, n_lo );
    more2 = false;

    compnz_next ( n2, k, a, more2 );
  }
//
//  Do the next element at the current value of N.
//
  else if ( more2 )
  {
    compnz_next ( n2, k, a, more2 );
  }
  else
  {
    n2 = n2 + 1;

    compnz_next ( n2, k, a, more2 );
  }
//
//  Termination occurs if MORE2 = FALSE and N2 = N_HI.
//
  if ( !more2 && n2 == n_hi )
  {
    more = false;
  }
  return;
}
//****************************************************************************80

void subset_by_size_next ( int n, int a[], int &size, bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_BY_SIZE_NEXT returns all subsets of an N set, in order of size.
//
//  Example:
//
//    N = 4:
//
//    1 2 3 4
//    1 2 3
//    1 2 4
//    1 3 4
//    1 3
//    1 4
//    2 3
//    1
//    2
//    3
//    (the empty set)
//
//  Discussion:
//
//    The subsets are returned in decreasing order of size, with the
//    empty set last.
//
//    For a given size K, the K subsets are returned in lexicographic order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the set.
//
//    Input/output, int A[N].  The entries A(1:SIZE) contain
//    the elements of the subset.  The elements are given in ascending
//    order.
//
//    Output, int &SIZE, the number of elements in the subset.
//
//    Input/output, bool &MORE.  Set MORE = FALSE before first call
//    for a new sequence of subsets.  It then is set and remains
//    TRUE as long as the subset computed on this call is not the
//    final one.  When the final subset is computed, MORE is set to
//    FALSE as a signal that the computation is done.
//
{
  static bool more2 = false;

  if ( !more )
  {
    more = true;
    more2 = false;
    size = n;
  }
  else if ( !more2 )
  {
    size = size - 1;
  }
//
//  Compute the next subset of size SIZE.
//
  if ( 0 < size )
  {
    ksub_next ( n, size, a, more2 );
  }
  else if ( size == 0 )
  {
    more = false;
  }
  return;
}
//****************************************************************************80

void subset_gray_next ( int n, int a[], bool &more, int &ncard, int &iadd )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_GRAY_NEXT generates all subsets of a set of order N, one at a time.
//
//  Discussion:
//
//    It generates the subsets one at a time, by adding or subtracting
//    exactly one element on each step.
//
//    The user should set MORE = .FALSE. and the value of N before
//    the first call.  On return, the user may examine A which contains
//    the definition of the new subset, and must check .MORE., because
//    as soon as it is .FALSE. on return, all the subsets have been
//    generated and the user probably should cease calling.
//
//    The first set returned is the empty set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the order of the total set from which
//    subsets will be drawn.
//
//    Input/output, int A[N].  On each return, the Gray code for the newly
//    generated subset.  A[I] = 0 if element I is in the subset, 1 otherwise.
//
//    Input/output, bool &MORE.  Set this variable FALSE before
//    the first call.  Normally, MORE will be returned TRUE but once
//    all the subsets have been generated, MORE will be
//    reset FALSE on return and you should stop calling the program.
//
//    Input/output, int &NCARD, the cardinality of the set returned,
//    which may be any value between 0 (the empty set) and N (the
//    whole set).
//
//    Output, int &IADD, the element which was added or removed to the
//    previous subset to generate the current one.  Exception:
//    the empty set is returned on the first call, and IADD is set to -1.
{
  int i;
//
//  First set returned is the empty set.
//
  if ( !more )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0;
    }

    iadd = 0;
    ncard = 0;
    more = true;
  }
  else
  {
    iadd = 1;

    if ( ( ncard % 2 ) != 0 )
    {
      for ( ; ; )
      {
        iadd = iadd + 1;
        if ( a[iadd-2] != 0 )
        {
          break;
        }
      }
    }

    a[iadd-1] = 1 - a[iadd-1];
    ncard = ncard + 2 * a[iadd-1] - 1;
//
//  Last set returned is the singleton A(N).
//
    if ( ncard == a[n-1] )
    {
      more = false;
    }
  }
  return;
}
//****************************************************************************80

int subset_gray_rank ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_GRAY_RANK ranks a subset of an N set, using the Gray code ordering.
//
//  Example:
//
//    N = 4
//
//       A       Rank
//    -------   -----
//
//    0 0 0 0       1
//    0 0 0 1       2
//    0 0 1 1       3
//    0 0 1 0       4
//    0 1 1 0       5
//    0 1 1 1       6
//    0 1 0 1       7
//    0 1 0 0       8
//    1 1 0 0       9
//    1 1 0 1      10
//    1 1 1 1      11
//    1 1 1 0      12
//    1 0 1 0      13
//    1 0 1 1      14
//    1 0 0 1      15
//    1 0 0 0      16
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the total set from which
//    subsets will be drawn.
//
//    Input, int A[N]; A(I) is 1 if element I is in the set,
//    and 0 otherwise.
//
//    Output, int SUBSET_GRAY_RANK, the rank of the subset in the 
//    Gray code ordering.
//
{
  int gray;
  int rank;

  gray = ( int ) ubvec_to_ui4 ( n, a );

  rank = gray_rank ( gray );

  rank = rank + 1;

  return rank;
}
//****************************************************************************80

void subset_gray_unrank ( int rank, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_GRAY_UNRANK produces a subset of an N set of the given Gray code rank.
//
//  Example:
//
//    N = 4
//
//     Rank     A    
//    -----  -------
//
//        1  0 0 0 0
//        2  0 0 0 1
//        3  0 0 1 1
//        4  0 0 1 0
//        5  0 1 1 0
//        6  0 1 1 1
//        7  0 1 0 1
//        8  0 1 0 0
//        9  1 1 0 0
//       10  1 1 0 1
//       11  1 1 1 1
//       12  1 1 1 0
//       13  1 0 1 0
//       14  1 0 1 1
//       15  1 0 0 1
//       16  1 0 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int RANK, the rank of the subset in the Gray code ordering.
//
//    Input, int N, the order of the total set from which
//    subsets will be drawn.
//
//    Output, int A[N]; A(I) is 1 if element I is in the set,
//    and 0 otherwise.
//
{
  int gray;

  gray = ( unsigned int ) gray_unrank ( rank-1 );

  ui4_to_ubvec ( gray, n, a );

  return;
}
//****************************************************************************80

void subset_lex_next ( int n, bool jmp, int ndim, int &k, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_LEX_NEXT generates the subsets of a set of N elements, one at a time.
//
//  Discussion:
//
//    The subsets are generated in lexicographical order.  
//
//    The routine can also be forced to generate only those subsets whose 
//    size is no greater than some user-specified maximum.
//
//  Example:
//
//    N = 5, JMP = ( K == 3 )
//
//    1
//    1 2
//    1 2 3
//    1 2 4
//    1 2 5
//    1 3
//    1 3 4
//    1 3 5
//    1 4
//    1 4 5
//    1 5
//    2
//    2 3
//    2 3 4
//    2 3 5
//    2 4
//    2 4 5
//    2 5
//    3
//    3 4
//    3 4 5
//    3 5
//    4
//    4 5
//    5
//    empty set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the order of the main set from which subsets
//    are chosen.
//
//    Input, bool JMP.  In the simplest case, set JMP = .FALSE. for
//    a normal computation.  But to jump over supersets of the input set,
//    set JMP = TRUE.  Setting JMP = ( K == 3 ) before every new call
//    will, for example, force all the subsets returned
//    to have cardinality 3 or less.
//
//    Input, int NDIM, the allowed storage for A.  If NDIM < N,
//    JMP must be used to avoid creation of a subset too large to store in A.
//
//    Input/output, int &K.  On first call, the user must set K = 0 as
//    a startup signal to the program.  Thereafter, the routine returns
//    the size of the computed subset in K.  On the last return,
//    the empty set is returned and K is 0, which is a signal to
//    the user that the computation is complete.
//
//    Input/output, int A[NDIM].  A(I) is the I-th element of the
//    subset, listed in increasing order, with 0's in entries
//    beyond entry K.
//
{
  int is;

  if ( k <= 0 )
  {
    if ( jmp )
    {
      return;
    }
    is = 0;
    k = 1;
    a[0] = 1;
  }
  else if ( a[k-1] != n )
  {
    is = a[k-1];

    if ( !jmp )
    {
      k = k + 1;
    }

    a[k-1] = is + 1;
  }
  else
  {
    k = k - 1;

    if ( k != 0 )
    {
      a[k-1] = a[k-1] + 1;
    }
  }
  return;
}
//****************************************************************************80

void subset_random ( int n, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_RANDOM selects a random subset of an N-set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the full set.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int A[N].  A vector to hold the information about
//    the set chosen.  On return, if A[I] = 1, then
//    I is in the random subset, otherwise, A[I] = 0
//    and I is not in the random subset.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = i4_uniform ( 0, 1, seed );
  }
  return;
}
//****************************************************************************80

void subtriangle_next ( int n, bool &more, int &i1, int &j1, int &i2, int &j2, 
  int &i3, int &j3 )

//****************************************************************************80
//
//  Purpose:
//
//    SUBTRIANGLE_NEXT computes the next subtriangle of a triangle.
//
//  Discussion:
//
//    The three sides of a triangle have been subdivided into N segments,
//    inducing a natural subdivision of the triangle into N*N subtriangles.
//    It is desired to consider each subtriangle, one at a time, in some
//    definite order.  This routine can produce information defining each 
//    of the subtriangles, one after another.
//
//    The subtriangles are described in terms of the integer coordinates 
//    (I,J) of their vertices.  These coordinates both range from 0 to N,
//    with the additional restriction that I + J <= N.
//
//    The vertices of each triangle are listed in counterclockwise order.
//
//  Example:
//
//    N = 4
//
//    4  *
//       |\
//       16\
//    3  *--*
//       |14|\
//       13\15\
//    2  *--*--*
//       |\9|11|\
//       |8\10\12\
//    1  *--*--*--*
//       |\2|\4|\6|\
//       |1\|3\|5\|7\
//   0   *--*--*--*--*
//
//       0  1  2  3  4
//
//    Rank  I1 J1  I2 J2  I3 J3
//    ----  -----  -----  ----- 
//       1   0  0   1  0   0  1
//       2   1  1   0  1   1  0
//       3   1  0   2  0   1  1
//       4   2  1   1  1   2  0
//       5   2  0   3  0   2  1
//       6   3  1   1  1   3  0
//       7   3  0   4  0   3  1
//       8   0  1   1  1   0  2
//       9   1  2   0  2   1  1
//      10   1  1   2  1   1  2
//      11   2  2   1  2   2  1
//      12   2  1   3  1   2  2
//      13   0  2   1  2   0  3
//      14   1  3   0  3   1  2
//      15   1  2   2  2   1  3
//      16   0  3   1  3   0  4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, indicates the number of subdivisions of each side
//    of the original triangle.
//
//    Input/output, bool &MORE.
//    On first call, set MORE to FALSE.  Thereafter, the output value of MORE
//    will be TRUE if there are more subtriangles that can be generated by
//    further calls.  However, if MORE is returned as FALSE, the accompanying
//    subtriangle information refers to the last subtriangle that can be
//    generated.
//
//    Input/output, int &I1, &J1, &I2, &J2, &I3, &J3, the indices of the 
//    vertices of the subtriangle.
//
{
  if ( n <= 0 )
  {
    more = false;
    return;
  }

  if ( !more )
  {
    i1 = 0;
    j1 = 0;
    i2 = 1;
    j2 = 0;
    i3 = 0;
    j3 = 1;

    if ( n == 1 )
    {
      more = false;
    }
    else
    {
      more = true;
    }
  }
//
//  We last generated a triangle like:
//
//    2---1
//     \  |
//      \ |
//       \|
//        3
//
  else if ( i2 < i3 )
  {
    i1 = i3;
    j1 = j3;
    i2 = i1 + 1;
    j2 = j1;
    i3 = i1;
    j3 = j1 + 1;
  }
//
//  We last generated a triangle like
//
//    3
//    |\
//    | \
//    |  \
//    1---2
//
  else if ( i1 + 1 + j1 + 1 <= n )
  {
    i1 = i1 + 1;
    j1 = j1 + 1;
    i2 = i1 - 1;
    j2 = j1;
    i3 = i1;
    j3 = j1 - 1;
  }
//
//  We must be at the end of a row.
//
  else
  {
    i1 = 0;
    j1 = j1 + 1;
    i2 = i1 + 1;
    j2 = j1;
    i3 = i1;
    j3 = j1 + 1;

    if ( n <= j1 + 1 )
    {
      more = false;
    }
  }

  return;
}
//****************************************************************************80

void thue_binary_next ( int &n, int thue[] )

//****************************************************************************80
//
//  Purpose:
//
//    THUE_BINARY_NEXT returns the next element in a binary Thue sequence.
//
//  Discussion:
//
//    Thue demonstrated that arbitrarily long sequences of 0's and
//    1's could be generated which had the "cubefree" property.  In
//    other words, for a given string S, there was no substring W
//    such that S contained "WWW".  In fact, a stronger result holds:
//    if "a" is the first letter of W, it is never the case that S
//    contains the substring "WWa".
//
//    In this example, the digits allowed are binary, that is, just
//    "0" and "1".  The replacement rules are:
//
//    "0" -> "01"
//    "1" -> "10"
//
//    This routine produces the next binary Thue sequence in a given series.
//    However, the input sequence must be a Thue sequence in order for
//    us to guarantee that the output sequence will also have the
//    cubic nonrepetition property.
//
//    Also, enough space must be set aside in THUE to hold the
//    output sequence.  This will always be twice the input
//    value of N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &N.  On input, the length of the input sequence.
//    On output, the length of the output sequence.
//
//    Input, int THUE[N].  On input, the initial Thue sequence, and on
//    output, the result of applying the substitution rules once.
//
{
  int i;
  int n_out;
  int *thue_out;

  n_out = 0;
  thue_out = new int[2 * n ];

  for ( i = 0; i < n; i++ )
  {
    if ( thue[i] == 0 )
    {
      thue_out[n_out] = 0;
      n_out = n_out + 1;
      thue_out[n_out] = 1;
      n_out = n_out + 1;
    }
    else if ( thue[i] == 1 )
    {
      thue_out[n_out] = 1;
      n_out = n_out + 1;
      thue_out[n_out] = 0;
      n_out = n_out + 1;
    }
    else
    {
      cerr << "\n";
      cerr << "THUE_BINARY_NEXT - Fatal error!\n";
      cerr << "  The input sequence contains a non-binary digit\n";
      cerr << "  THUE[" << i << "] = " << thue[i] << "\n";
      exit ( 1 );
    }
  }

  n = n_out;

  for ( i = 0; i < n; i++ )
  {
    thue[i] = thue_out[i];
  }

  delete [] thue_out;

  return;
}
//****************************************************************************80

void thue_ternary_next ( int &n, int thue[] )

//****************************************************************************80
//
//  Purpose:
//
//    THUE_TERNARY_NEXT returns the next element in a ternary Thue sequence.
//
//  Discussion:
//
//    Thue was interested in showing that there were arbitrarily long
//    sequences of digits which never displayed a pair of contiguous
//    repetitions of any length.  That is, there was no occurrence of
//    "00" or "1010" or "121121", anywhere in the string.  This makes
//    the string "squarefree".
//
//    To do this, he demonstrated a way to start with a single digit,
//    and to repeatedly apply a series of transformation rules to each
//    digit of the sequence, deriving nonrepeating sequences of ever
//    greater length.
//
//    In this example, the digits allowed are ternary, that is, just
//    "0", "1" and "2".  The replacement rules are:
//
//    "0" -> "12"
//    "1" -> "102"
//    "2" -> "0"
//
//    This routine produces the next Thue sequence in a given series.
//    However, the input sequence must be a Thue sequence in order for
//    us to guarantee that the output sequence will also have the
//    nonrepetition property.
//
//    Also, enough space must be set aside in THUE to hold the
//    output sequence.  This will never be more than 3 times the input
//    value of N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Brian Hayes,
//    Third Base,
//    American Scientist,
//    Volume 89, Number 6, pages 490-494, November-December 2001.
//
//  Parameters:
//
//    Input/output, int &N.  On input, the length of the input sequence.
//    On output, the length of the output sequence.
//
//    Input, int THUE[*N].  On input, the initial Thue sequence, and on
//    output, the result of applying the substitution rules once.
//
{
  int i;
  int n_out;
  int *thue_out;

  n_out = 0;
  thue_out = new int[ 3 * n ];

  for ( i = 0; i < n; i++ )
  {

    if ( thue[i] == 0 )
    {
      thue_out[n_out] = 1;
      n_out = n_out + 1;
      thue_out[n_out] = 2;
      n_out = n_out + 1;
    }
    else if ( thue[i] == 1 )
    {
      thue_out[n_out] = 1;
      n_out = n_out + 1;
      thue_out[n_out] = 0;
      n_out = n_out + 1;
      thue_out[n_out] = 2;
      n_out = n_out + 1;
    }
    else if ( thue[i] == 2 )
    {
      thue_out[n_out] = 0;
      n_out = n_out + 1;
    }
    else
    {
      cerr << "\n";
      cerr << "THUE_TERNARY_NEXT - Fatal error!\n";
      cerr << "  The input sequence contains a non-ternary digit\n";
      cerr << "  THUE[" << i << "] = " << thue[i] << "\n";
      exit ( 1 );
    }
  }

  n = n_out;
  for ( i = 0; i < n_out; i++ )
  {
    thue[i] = thue_out[i];
  }

  delete [] thue_out;

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
//****************************************************************************80

void triang ( int n, int zeta[], int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANG renumbers elements in accordance with a partial ordering.
//
//  Discussion:
//
//    TRIANG is given a partially ordered set.  The partial ordering
//    is defined by a matrix ZETA, where element I is partially less than
//    or equal to element J if and only if ZETA(I,J) = 1.
//
//    TRIANG renumbers the elements with a permutation P so that if
//    element I is partially less than element J in the partial ordering,
//    then P(I) < P(J) in the usual, numerical ordering.
//
//    In other words, the elements are relabeled so that their labels
//    reflect their ordering.  This is equivalent to relabeling the
//    matrix so that, on unscrambling it, the matrix would be upper
//    triangular.
//
//    Calling I4MAT_PERM or R8MAT_PERM with P used for both the row
//    and column permutations applied to matrix ZETA will result in
//    an upper triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the number of elements in the set.
//
//    Input, int ZETA[N*N], describes the partial ordering.
//    ZETA[I,J] =:
//      0, for diagonal elements (I = J), or
//         for unrelated elements, or
//         if J << I.
//      1, if I << J.
//
//    Output, int P[N], a permutation of the elements that reflects
//    their partial ordering.  P[I] is the new label of element I, with
//    the property that if ZETA[I,J] = 1, that is, I << J,
//    then P[I] < P[J] (in the usual ordering).
//
{
  int i;
  bool error;
  int iq;
  int ir;
  int it;
  int l;
  int m;
//
//  Make sure ZETA represents a partially ordered set.  In other words,
//  if ZETA(I,J) = 1, then ZETA(J,I) must NOT be 1.
//
  error = pord_check ( n, zeta );

  if ( error )
  {
    cerr << "\n";
    cerr << "TRIANG - Fatal error!\n";
    cerr << "  The matrix ZETA does not represent a\n";
    cerr << "  partial ordering.\n";
    exit ( 1 );
  }

  m = 0;
  l = 0;
  for ( i = 0; i < n; i++ )
  {
    p[i] = 0;
  }
//
//  Find the next value of M for which P(M) is 0.
//
  for ( ; ; )
  {
    m = m + 1;

    if ( p[m-1] == 0 ) 
    {
      break;
    }

    if ( m == n )
    {
      return;
    }

  }

  it = m + 1;
  ir = m + 1;

  for ( ; ; )
  {
    if ( ir <= n )
    {
      if ( p[ir-1] == 0 && zeta[(ir-1)+(m-1)*n] != 0 )
      {
        p[ir-1] = m;
        m = ir;
        ir = it;
      }
      else
      {
        ir = ir + 1;
      }
    }
    else
    {
      l = l + 1;
      iq = p[m-1];
      p[m-1] = l;

      if ( iq != 0 )
      {

        ir = m + 1;
        m = iq;
      }
      else if ( m == n )
      {
        break;
      }
      else
      {
        for ( ; ; )
        {
          m = m + 1;

          if ( p[m-1] == 0 )
          {
            break;
          }

          if ( m == n )
          {
            return;
          }
        }
        it = m + 1;
        ir = m + 1;
      }
    }
  }
  return;
}
//****************************************************************************80

void tuple_next ( int m1, int m2, int n, int &rank, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT computes the next element of a tuple space.
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between M1 and M2.  The elements are produced one at a time.
//    The first element is
//      (M1,M1,...,M1),
//    the second element is
//      (M1,M1,...,M1+1),
//    and the last element is
//      (M2,M2,...,M2)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2, M1 = 1, M2 = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank  X      Rank   X
//    ----  ---    -----  ---
//    0     * *    1      1 1
//    1     1 1    2      1 2
//    2     1 2    3      1 3
//    3     1 3    4      2 1
//    4     2 1    5      2 2
//    5     2 2    6      2 3
//    6     2 3    7      3 1
//    7     3 1    8      3 2
//    8     3 2    9      3 3
//    9     3 3    0      0 0
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
//    Input, int M1, M2, the minimum and maximum entries.
//
//    Input, int N, the number of components.
//
//    Input/output, int &RANK, counts the elements.
//    On first call, set RANK to 0.  Thereafter, the output value of RANK
//    will indicate the order of the element returned.  When there are no
//    more elements, RANK will be returned as 0.
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
{
  int i;
  int j;

  if ( m2 < m1 )
  {
    rank = 0;
    return;
  }

  if ( rank <= 0 )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = m1;
    }
    rank = 1;
  }
  else
  {
    rank = rank + 1;
    i = n - 1;

    for ( ; ; )
    {

      if ( x[i] < m2 )
      {
        x[i] = x[i] + 1;
        break;
      }

      x[i] = m1;

      if ( i == 0 )
      {
        rank = 0;
        for ( j = 0; j < n; j++ )
        {
          x[j] = m1;
        }
        break;
      }
      i = i - 1;
    }
  }
  return;
}
//****************************************************************************80

void tuple_next_fast ( int m, int n, int rank, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between 1 and M.  The elements are produced one at a time.
//    The first element is
//      (1,1,...,1)
//    and the last element is
//      (M,M,...,M)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2,
//    M = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank          X
//    ----          ----
//   -1            -1 -1
//
//    0             1  1
//    1             1  2
//    2             1  3
//    3             2  1
//    4             2  2
//    5             2  3
//    6             3  1
//    7             3  2
//    8             3  3
//    9             1  1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the maximum entry in each component.
//    M must be greater than 0.
//
//    Input, int N, the number of components.
//    N must be greater than 0.
//
//    Input, integer RANK, indicates the rank of the tuples.
//    Typically, 0 <= RANK < N**M; values larger than this are legal
//    and meaningful, and are equivalent to the corresponding value
//    MOD N**M.  If RANK < 0, this indicates that this is the first
//    call for the given values of (M,N).  Initialization is done,
//    and X is set to a dummy value.
//
//    Output, int X[N], the next tuple, or a dummy value if initialization
//    is being done.
//
{
  static int *base = NULL;
  int i;

  if ( rank < 0 )
  {
    if ( m <= 0 )
    {
      cerr << "\n";
      cerr << "TUPLE_NEXT_FAST - Fatal error!\n";
      cerr << "  The value M <= 0 is not legal.\n";
      cerr << "  M = " << m << "\n";
      exit ( 1 );
    }
    if ( n <= 0 )
    {
      cerr << "\n";
      cerr << "TUPLE_NEXT_FAST - Fatal error!\n";
      cerr << "  The value N <= 0 is not legal.\n";
      cerr << "  N = " << n << "\n";
      exit ( 1 );
    }

    if ( base )
    {
      delete [] base;
    }
    base = new int[n];

    base[n-1] = 1;
    for ( i = n-2; 0 <= i; i-- )
    {
      base[i] = base[i+1] * m;
    }
    for ( i = 0; i < n; i++ )
    {
      x[i] = -1;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( rank / base[i] ) % m ) + 1;
    }
  }
  return;
}
//****************************************************************************80

void tuple_next_ge ( int m, int n, int &k, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT_GE computes the next "nondecreasing" element of a tuple space.
//
//  Discussion:
//
//    The elements are N vectors.  Each element is constrained to lie
//    between 1 and M, and to have components that are nondecreasing.
//    That is, for an element X, and any positive K,
//      X(I) <= X(I+K)
//
//    The elements are produced one at a time.
//    The first element is
//      (1,1,...,1)
//    and the last element is
//      (M,M,...,M)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 3, M = 3
//
//    K   X
//    --  -----
//     1  1 1 1
//     2  1 1 2
//     3  1 1 3
//     4  1 2 2
//     5  1 2 3
//     6  1 3 3
//     7  2 2 2
//     8  2 2 3
//     9  2 3 3
//    10  3 3 3
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
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int M, the maximum entry.
//
//    Input, int N, the number of components.
//
//    Input/output, int &K, counts the elements.
//    On first call, set K to 0.  Thereafter, K will indicate the
//    order of the element returned.  When there are no more elements,
//    K will be returned as 0.
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
{
  int i;
  int j;

  if ( m < 1 )
  {
    return;
  }

  if ( k <= 0 )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = 1;
    }
    k = 1;
    return;
  }

  for ( i = n-1; 0 <= i; i-- )
  {
    if ( x[i] < m )
    {
      x[i] = x[i] + 1;
      for ( j = i+1; j < n; j++ )
      {
        x[j] = x[i];
      }
      k = k + 1;
      return;
    }
  }

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    x[j] = 0;
  }

  return;
}
//****************************************************************************80

void tuple_next2 ( int n, int xmin[], int xmax[], int x[], int &rank )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT2 computes the next element of an integer tuple space.
//
//  Discussion:
//
//    The elements X are N vectors.
//
//    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
//
//    The elements are produced one at a time.
//
//    The first element is
//      (XMIN(1), XMIN(2), ..., XMIN(N)),
//    the second is (probably)
//      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
//    and the last element is
//      (XMAX(1), XMAX(2), ..., XMAX(N))
//
//    Intermediate elements are produced in a lexicographic order, with
//    the first index more important than the last, and the ordering of
//    values at a fixed index implicitly defined by the sign of
//    XMAX(I) - XMIN(I).
//
//  Example:
//
//    N = 2,
//    XMIN = (/ 1, 10 /)
//    XMAX = (/ 3,  8 /)
//
//    RANK    X
//    ----  -----
//      1   1 10
//      2   1  9
//      3   1  8
//      4   2 10
//      5   2  9
//      6   2  8
//      7   3 10
//      8   3  9
//      9   3  8
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
//    Input, int N, the number of components.
//
//    Input, int XMIN[N], XMAX[N], the "minimum" and "maximum" entry values.
//    These values are minimum and maximum only in the sense of the lexicographic
//    ordering.  In fact, XMIN(I) may be less than, equal to, or greater
//    than XMAX(I).
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
//    Input/output, int &RANK, the rank of the item.  On first call,
//    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
//    there are no more items in the sequence.
//
{
  int i;
  int test;

  if ( rank < 0 )
  {
    cerr << "\n";
    cerr << "TUPLE_NEXT2 - Fatal error!\n";
    cerr << "  Illegal value of RANK = " << rank << "\n";
    exit ( 1 );
  }

  test = 1;
  for ( i = 0; i < n; i++ )
  {
    test = test * ( 1 + abs ( xmax[i] - xmin[i] ) );
  }

  if ( test < rank )
  {
    cerr << "\n";
    cerr << "TUPLE_NEXT2 - Fatal error!\n";
    cerr << "  Illegal value of RANK = " << rank << "\n";
    exit ( 1 );
  }

  if ( rank == 0 )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = xmin[i];
    }
    rank = 1;
    return;
  }

  rank = rank + 1;
  i = n - 1;

  for ( ; ; )
  {
    if ( x[i] != xmax[i] )
    {
      if ( xmin[i] < xmax[i] )
      {
        x[i] = x[i] + 1;
      }
      else
      {
        x[i] = x[i] - 1;
      }
      break;
    }

    x[i] = xmin[i];

    if ( i == 0 )
    {
      rank = 0;
      break;
    }

    i = i - 1;
  }

  return;
}
//****************************************************************************80

bool ubvec_add ( int n, int bvec1[], int bvec2[], int bvec3[] )

//****************************************************************************80
//
//  Purpose:
//
//    UBVEC_ADD adds two unsigned binary vectors.
//
//  Discussion:
//
//    A UBVEC is a vector of binary digits representing an unsigned integer.  
//
//    UBVEC[N-1] contains the units digit, UBVEC[N-2]
//    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Example:
//
//    N = 4
//
//     UBVEC1       +  UBVEC2       =  UBVEC3
//
//    ( 0 0 0 1 )   + ( 0 0 1 1 )   = ( 0 1 0 0 )
//
//      1           +   3           =   4
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
//    Output, bool UBVEC_ADD, is TRUE if an error occurred.
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

unsigned int ubvec_to_ui4 ( int n, int bvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    UBVEC_TO_UI4 makes an unsigned integer from an unsigned binary vector.
//
//  Discussion:
//
//    A UBVEC is a vector of binary digits representing an unsigned integer.  
//
//    UBVEC[N-1] contains the units digit, UBVEC[N-2]
//    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//  Example:
//
//    N = 4
//
//         BVEC   binary  I
//    ----------  -----  --
//    1  2  3  4
//    ----------
//    0  0  0  1       1  1
//    0  0  1  0      10  2
//    0  0  1  1      11  3
//    0  1  0  0     100  4
//    1  0  0  1    1001  9
//    1  1  1  1    1111 15
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
//    Input, unsigned int UBVEC_TO_UI4, the integer value.
//
{
  int base = 2;
  int i;
  unsigned int ui4;

  ui4 = 0;
  for ( i = 0; i < n; i++ )
  {
    ui4 = base * ui4 + bvec[i];
  }

  return ui4;
}
//****************************************************************************80

void ui4_to_ubvec ( unsigned int ui4, int n, int bvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    UI4_TO_UBVEC makes a unsigned binary vector from an integer.
//
//  Discussion:
//
//    A UBVEC is a vector of binary digits representing an unsigned integer.  
//
//    UBVEC[N-1] contains the units digit, UBVEC[N-2]
//    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
//    so that printing the digits in order gives the binary form of the number.
//
//    To guarantee that there will be enough space for any
//    value of I, it would be necessary to set N = 32.
//
//  Example:
//
//    UI4      BVEC         binary
//        0  1  2  3  4  5
//        1  2  4  8 16 32
//    --  ----------------  ------
//     1  1  0  0  0  0  1       1
//     2  0  1  0  0  1  0      10
//     3  1  1  0  0  1  1      11
//     4  0  0  0  1  0  0     100
//     9  0  0  1  0  0  1    1001
//    57  1  1  1  0  1  1  110111
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
//    Input, unsigned int UI4, an integer to be represented.
//
//    Input, int N, the dimension of the vector.
//
//    Output, int BVEC[N], the binary representation.
//
{
  int base = 2;
  int i;

  for ( i = n - 1; 0 <= i; i-- )
  {
    bvec[i] = ui4 % base;
    ui4 = ui4 / base;
  }

  return;
}
//****************************************************************************80

void vec_colex_next ( int dim_num, int base, int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_COLEX_NEXT generates vectors in colex order.
//
//  Discussion:
//
//    The vectors are produced in colexical order, starting with
//    (0,0,...,0),
//    (1,0,...,0), 
//    ...
//    (BASE-1,BASE-1,...,BASE-1).
//
//  Example:
//
//    DIM_NUM = 2,
//    BASE = 3
//
//    0   0
//    1   0
//    2   0
//    0   1
//    1   1
//    2   1
//    0   2
//    1   2
//    2   2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int BASE, the base to be used.  BASE = 2 will
//    give vectors of 0's and 1's, for instance.
//
//    Output, int A[DIM_NUM], the next vector.
//
//    Input/output, bool &MORE.  Set this variable false before
//    the first call.  On return, MORE is TRUE if another vector has
//    been computed.  If MORE is returned FALSE, ignore the output 
//    vector and stop calling the routine.
//
{
  int i;

  if ( !more )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = 0;
    }
    more = true;
  }
  else
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = a[i] + 1;

      if ( a[i] < base )
      {
        return;
      }
      a[i] = 0;
    }
    more = false;
  }

  return;
}
//****************************************************************************80

void vec_colex_next2 ( int dim_num, int base[], int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_COLEX_NEXT2 generates vectors in colex order.
//
//  Discussion:
//
//    The vectors are produced in colexical order, starting with
//
//    (0,        0,        ...,0),
//    (1,        0,        ...,0),
//     ...
//    (BASE(1)-1,0,        ...,0)
//
//    (0,        1,        ...,0)
//    (1,        1,        ...,0)
//    ...
//    (BASE(1)-1,1,        ...,0)
//
//    (0,        2,        ...,0)
//    (1,        2,        ...,0)
//    ...
//    (BASE(1)-1,BASE(2)-1,...,BASE(DIM_NUM)-1).
//
//  Example:
//
//    DIM_NUM = 2,
//    BASE = { 3, 3 }
//
//    0   0
//    1   0
//    2   0
//    0   1
//    1   1
//    2   1
//    0   2
//    1   2
//    2   2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int BASE[DIM_NUM], the bases to be used in each dimension.
//    In dimension I, entries will range from 0 to BASE[I]-1.
//
//    Output, int A[DIM_NUM], the next vector.
//
//    Input/output, bool &MORE.  Set this variable false before
//    the first call.  On return, MORE is TRUE if another vector has
//    been computed.  If MORE is returned FALSE, ignore the output 
//    vector and stop calling the routine.
//
{
  int i;

  if ( !more )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = 0;
    }
    more = true;
  }
  else
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = a[i] + 1;

      if ( a[i] < base[i] )
      {
        return;
      }
      a[i] = 0;
    }
    more = false;
  }

  return;
}
//****************************************************************************80

void vec_colex_next3 ( int dim_num, int base[], int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_COLEX_NEXT3 generates vectors in colex order.
//
//  Discussion:
//
//    The vectors are produced in colexical order, starting with
//
//    (1,        1,        ...,1),
//    (2,        1,        ...,1),
//     ...
//    (BASE(1),  1,        ...,1)
//
//    (1,        2,        ...,1)
//    (2,        2,        ...,1)
//    ...
//    (BASE(1),  2,        ...,1)
//
//    (1,        3,        ...,1)
//    (2,        3,        ...,1)
//    ...
//    (BASE(1),  BASE(2), ...,BASE(DIM_NUM)).
//
//  Example:
//
//    DIM_NUM = 2,
//    BASE = { 3, 3 }
//
//    1   1
//    2   1
//    3   1
//    1   2
//    2   2
//    3   2
//    1   3
//    2   3
//    3   3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int BASE[DIM_NUM], the bases to be used in each dimension.
//    In dimension I, entries will range from 1 to BASE[I].
//
//    Output, int A[DIM_NUM], the next vector.
//
//    Input/output, bool &MORE.  Set this variable false before
//    the first call.  On return, MORE is TRUE if another vector has
//    been computed.  If MORE is returned FALSE, ignore the output 
//    vector and stop calling the routine.
//
{
  int i;

  if ( !more )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = 1;
    }
    more = true;
  }
  else
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = a[i] + 1;

      if ( a[i] <= base[i] )
      {
        return;
      }
      a[i] = 1;
    }
    more = false;
  }
  return;
}
//****************************************************************************80

void vec_gray_next ( int n, int base[], int a[], bool &done, int &change )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_GRAY_NEXT computes the elements of a product space.
//
//  Discussion:
//
//    The elements are produced one at a time.
//
//    This routine handles the case where the number of degrees of freedom may
//    differ from one component to the next.
//
//    A method similar to the Gray code is used, so that successive
//    elements returned by this routine differ by only a single element.
//
//    The routine uses internal static memory.
//
//  Example:
//
//    N = 2, BASE = ( 2, 3 ), DONE = TRUE
//
//     A    DONE  CHANGE
//    ---  -----  ------
//    0 0  FALSE    0 (1)  <-- C++ routine returns 0-based CHANGE.
//    0 1  FALSE    1 (2)
//    0 2  FALSE    1 (2)
//    1 2  FALSE    0 (1)
//    1 1  FALSE    1 (2)
//    1 0  FALSE    1 (2)
//    1 0   TRUE   -1  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, int BASE[N], contains the number of degrees of
//    freedom of each component.  The output values of A will
//    satisfy 0 <= A[I] < BASE[I].
//
//    Input/output, int A[N].  On the first call, the input value
//    of A doesn't matter.  Thereafter, it should be the same as
//    its output value from the previous call.  On output, if DONE
//    is FALSE, then A contains the next element of the space.
//
//    Input/output, bool &DONE.  On the first call, the user must
//    set DONE to TRUE.  This signals the program to initialize data.
//    On every return, if DONE is FALSE, the program has computed
//    another entry, which is contained in A.  If DONE is TRUE,
//    then there are no more entries, and the program should not be
//    called for any more.
//
//    Output, int &CHANGE, is set to the index of the element whose
//    value was changed.  On return from the first call, CHANGE
//    is 0, even though all the elements have been "changed".  On
//    return with DONE equal to TRUE, CHANGE is -1.  (Note that CHANGE
//    is a vector index.  In this C++ version, it is zero-based.)
//
{
  static int *active = NULL;
  int i;
  static int *dir = NULL;
//
//  The user is calling for the first time.
//
  if ( done )
  {
    done = false;
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0;
    }

    if ( active )
    {
      delete [] active;
    }
    active = new int[n];

    if ( dir )
    {
      delete [] dir;
    }
    dir = new int[n];

    for ( i = 0; i < n; i++ )
    {
      dir[i] = 1;
    }
    for ( i = 0; i < n; i++ )
    {
      active[i] = 1;
    }

    for ( i = 0; i < n; i++ )
    {
      if ( base[i] < 1 )
      {
        cerr << "\n";
        cerr << "VEC_GRAY_NEXT - Warning\n";
        cerr << "  For index I = " << i << "\n";
        cerr << "  the nonpositive value of BASE[I] = " << base[i] << "\n";
        cerr << "  which was reset to 1!\n";
        base[i] = 1;
        active[i] = 0;
      }
      else if ( base[i] == 1 )
      {
        active[i] = 0;
      }
    }

    change = 0;
    return;
  }
//
//  Find the maximum active index.
//
  change = -1;

  for ( i = 0; i < n; i++ )
  {
    if ( active[i] != 0 )
    {
      change = i;
    }
  }
//
//  If there are NO active indices, we have generated all vectors.
//
  if ( change == -1 )
  {
    delete [] dir;
    dir = NULL;
    delete [] active;
    active = NULL;
    done = true;
    return;
  }
//
//  Increment the element with maximum active index.
//
  a[change] = a[change] + dir[change];
//
//  If we attained a minimum or maximum value, reverse the direction
//  vector, and deactivate the index.
//
  if ( a[change] == 0 || a[change] == base[change] - 1 )
  {
    dir[change] = -dir[change];
    active[change] = 0;
  }
//
//  Activate all subsequent indices.
//
  for ( i = change + 1; i < n; i++ )
  {
    if ( 1 < base[i] )
    {
      active[i] = 1;
    }
  }

  return;
}
//****************************************************************************80

int vec_gray_rank ( int n, int base[], int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_GRAY_RANK computes the rank of a product space element.
//
//  Discussion:
//
//    The rank applies only to the elements as produced by the routine
//    VEC_GRAY_NEXT.
//
//  Example:
//
//    N = 2, BASE = ( 2, 3 ), A = ( 1, 2 ),
//
//    RANK = 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, int BASE[N], contains the number of degrees of
//    freedom of each component.  The output values of A will
//    satisfy 0 <= A[I] < BASE[I].
//
//    Input, int A[N], the product space element, with the
//    property that 0 <= A[I] < BASE[I] for each entry I.
//
//    Output, int VEC_RANK, the rank, or order, of the element in
//    the list of all elements.  The rank count begins at 1.
//
{
  int c;
  int i;
  int rank;

  rank = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( ( rank % 2 ) == 1 )
    {
      c = base[i] - a[i] - 1;
    }
    else
    {
      c = a[i];
    }
    rank = base[i] * rank + c;
  }

  rank = rank + 1;

  return rank;
}
//****************************************************************************80

void vec_gray_unrank ( int n, int base[], int rank, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_GRAY_UNRANK computes the product space element of a given rank.
//
//  Discussion:
//
//    The rank applies only to the elements as produced by the routine
//    VEC_GRAY_NEXT.
//
//  Example:
//
//    N = 2, BASE = ( 2, 3 ), RANK = 4.
//
//    A = ( 1, 2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, int BASE[N], contains the number of degrees of
//    freedom of each component.  The output values of A will
//    satisfy 0 <= A[I] < BASE[I].
//
//    Input, int RANK, the desired rank, or order, of the element in
//    the list of all elements.  The rank count begins at 1 and extends
//    to MAXRANK = Product ( 0 <= I <= N ) BASE[I].
//
//    Output, int A[N], the product space element of the given rank.
//
{
  int i;
  int s;

  s = rank - 1;

  for ( i = n-1; 0 <= i; i-- )
  {
    a[i] = s % base[i];
    s = s / base[i];

    if ( ( s % 2 ) == 1 )
    {
      a[i] = base[i] - a[i] - 1;
    }
  }

  return;
}
//****************************************************************************80

void vec_lex_next ( int dim_num, int base, int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_LEX_NEXT generates vectors in lex order.
//
//  Discussion:
//
//    The vectors are produced in lexical order, starting with
//    (0,0,...,0),
//    (0,0,...,1), 
//    ...
//    (BASE-1,BASE-1,...,BASE-1).
//
//  Example:
//
//    DIM_NUM = 2,
//    BASE = 3
//
//    0   0
//    0   1
//    0   2
//    1   0
//    1   1
//    1   2
//    2   0
//    2   1
//    2   2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the size of the vectors to be used.
//
//    Input, int BASE, the base to be used.  BASE = 2 will
//    give vectors of 0's and 1's, for instance.
//
//    Output, int A[DIM_NUM], the next vector.  
//
//    Input/output, bool &MORE.  Set this variable false before
//    the first call.  On return, MORE is TRUE if another vector has
//    been computed.  If MORE is returned FALSE, ignore the output 
//    vector and stop calling the routine.
//
{
  int i;

  if ( !more )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = 0;
    }
    more = true;
  }
  else
  {
    for ( i = dim_num - 1; 0 <= i; i-- )
    {
      a[i] = a[i] + 1;

      if ( a[i] < base )
      {
        return;
      }
      a[i] = 0;
    }
    more = false;
  }

  return;
}
//****************************************************************************80

void vec_random ( int n, int base, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_RANDOM selects a random N-vector of integers modulo a given base.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the vector to be generated.
//
//    Input, int BASE, the base to be used.
//
//    Input/output, int &SEED, a random number seed.
//
//    Output, int A[N], a list of N random values between
//    0 and BASE-1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = i4_uniform ( 0, base-1, seed );
  }

  return;
}
//****************************************************************************80

void vector_constrained_next ( int n, int x_min[], int x_max[], int x[], 
  int &constraint, bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VECTOR_CONSTRAINED_NEXT returns the "next" constrained vector.
//
//  Discussion:
//
//    We consider all vectors of dimension N whose components
//    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
//
//    We are only interested in the subset of these vectors which
//    satisfy the following constraint:
//
//      sum ( 1 <= I <= N ) ( ( X(I) - 1 ) / X_MAX(I) ) <= 1
//
//    We can carry out this check using integer arithmetic if we
//    multiply through by P = product ( X_MAX(1:N) ):
//
//      sum ( 1 <= I <= N ) ( ( X(I) - 1 ) * ( P / X_MAX(I) ) ) <= P.
//
//    This routine returns, one at a time, and in right-to-left
//    lexicographic order, exactly those vectors which satisfy
//    the constraint.
//
//  Example:
//
//    N = 3
//    X_MIN:   2   2   1
//    X_MAX:   4   5   3
//
//    P = 60
//
//    #  X(1)  X(2)  X(3)  CONSTRAINT
//
//    1    2     2     1       27
//    2    3     2     1       42
//    3    4     2     1       57
//    4    2     3     1       39
//    5    3     3     1       54
//    6    2     4     1       51
//    7    2     2     2       47
//    8    2     3     2       59
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components in the vector.
//
//    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
//    values allowed in each component.
//
//    Input/output, integer X[N].  On first call (with MORE = FALSE),
//    the input value of X is not important.  On subsequent calls, the
//    input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Output, int &CONSTRAINT, the constraint value for X.  Valid vectors X
//    will have a CONSTRAINT value between product(X_MIN(1:N)) (automatically)
//    and product(X_MAX(1:N)) (because we skip over vectors with a
//    constraint larger than this value).
//
//    Input/output, bool &MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  int i;
  int j;
  static int x_prod = -1;

  if ( !more )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = x_min[j];
    }

    x_prod = 1;
    for ( j = 0; j < n; j++ )
    {
      x_prod = x_prod * x_max[j];
    }

    constraint = 0;
    for ( j = 0; j < n; j++ )
    {
      constraint = constraint + ( ( x[j] - 1 ) * ( x_prod / x_max[j] ) );
    }

    if ( x_prod < constraint )
    {
      more = false;
    }
    else
    {
      more = true;
    }

    return;
  }
  else
  {
    i = 0;

    for ( ; ; )
    {
      if ( x[i] < x_max[i] )
      {
        x[i] = x[i] + 1;

        constraint = 0;
        for ( j = 0; j < n; j++ )
        {
          constraint = constraint + ( ( x[j] - 1 ) * ( x_prod / x_max[j] ) );
        }

        if ( constraint <= x_prod )
        {
          break;
        }
      }

      x[i] = x_min[i];

      i = i + 1;

      if ( n <= i )
      {
        more = false;
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

void vector_constrained_next2 ( int n, int x_min[], int x_max[], int x[], 
  int &constraint, bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VECTOR_CONSTRAINED_NEXT2 returns the "next" constrained vector.
//
//  Discussion:
//
//    We consider all vectors of dimension N whose components
//    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
//
//    We are only interested in the subset of these vectors which
//    satisfy the following constraint:
//
//      sum ( 1 <= I <= N ) ( X(I) / X_MAX(I) ) <= 1
//
//    We can carry out this check using integer arithmetic if we
//    multiply through by P = product ( X_MAX(1:N) ):
//
//      sum ( 1 <= I <= N ) ( X(I)  * ( P / X_MAX(I) ) ) <= P.
//
//    This routine returns, one at a time, and in right-to-left
//    lexicographic order, exactly those vectors which satisfy
//    the constraint.
//
//  Example:
//
//    N = 3
//    X_MIN:   1   1   1
//    X_MAX:   5   6   4
//
//    P = 120
//
//    #  X(1)  X(2)  X(3)  CONSTRAINT
//
//    1    1     1     1       74
//    2    2     1     1       98
//    3    1     2     1       94
//    4    2     2     1      119
//    5    1     3     1      114
//    6    1     1     2      104
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components in the vector.
//
//    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
//    values allowed in each component.
//
//    Input/output, integer X[N].  On first call (with MORE = FALSE),
//    the input value of X is not important.  On subsequent calls, the
//    input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Output, int &CONSTRAINT, the constraint value for X.  Valid vectors X
//    will have a CONSTRAINT value between product(X_MIN(1:N)) (automatically)
//    and product(X_MAX(1:N)) (because we skip over vectors with a
//    constraint larger than this value).
//
//    Input/output, bool &MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  int i;
  int j;
  static int x_prod = -1;

  if ( !more )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = x_min[j];
    }

    x_prod = 1;
    for ( j = 0; j < n; j++ )
    {
      x_prod = x_prod * x_max[j];
    }

    constraint = 0;
    for ( j = 0; j < n; j++ )
    {
      constraint = constraint + ( x[j] * ( x_prod / x_max[j] ) );
    }

    if ( x_prod < constraint )
    {
      more = false;
    }
    else
    {
      more = true;
    }

    return;
  }
  else
  {
    i = 0;

    for ( ; ; )
    {
      if ( x[i] < x_max[i] )
      {
        x[i] = x[i] + 1;

        constraint = 0;
        for ( j = 0; j < n; j++ )
        {
          constraint = constraint + ( x[j] * ( x_prod / x_max[j] ) );
        }

        if ( constraint <= x_prod )
        {
          break;
        }
      }

      x[i] = x_min[i];

      i = i + 1;

      if ( n <= i )
      {
        more = false;
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

void vector_constrained_next3 ( int n, int x_min[], int x_max[], int x[], 
  double &constraint, bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VECTOR_CONSTRAINED_NEXT3 returns the "next" constrained vector.
//
//  Discussion:
//
//    This routine addresses the same problem as VECTOR_CONSTRAINED_NEXT2,
//    and differs only in that real arithmetic is used, rather than
//    integer arithmetic.  Integer arithmetic allows us to do an exact
//    calculation, but we run into overflow problems in simple cases
//    where N is 10 and the X_MAX entries are of order 10, for instance.
//
//    We consider all vectors of dimension N whose components
//    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
//
//    We are only interested in the subset of these vectors which
//    satisfy the following constraint:
//
//      sum ( 1 <= I <= N ) ( X(I) / X_MAX(I) ) <= 1
//
//  Example:
//
//    N = 3
//    X_MIN:   1   1   1
//    X_MAX:   5   6   4
//
//    P = 120
//
//    #  X(1)  X(2)  X(3)  CONSTRAINT
//
//    1    1     1     1       0.62
//    2    2     1     1       0.82
//    3    1     2     1       0.78
//    4    2     2     1       0.98
//    5    1     3     1       0.95
//    6    1     1     2       0.87
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components in the vector.
//
//    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
//    values allowed in each component.
//
//    Input/output, integer X[N].  On first call (with MORE = FALSE),
//    the input value of X is not important.  On subsequent calls, the
//    input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Output, double &CONSTRAINT, the constraint value for X.  Valid vectors 
//    X will have a CONSTRAINT value between 
//      product(X_MIN(1:N)) / product(X_MAX(1:N))
//    and 1.0.
//
//    Input/output, bool &MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  int i;
  int j;

  if ( !more )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = x_min[j];
    }

    constraint = 0.0;
    for ( j = 0; j < n; j++ )
    {
      constraint = constraint + ( double ) x[j] / ( double ) x_max[j];
    }

    if ( 1.0 < constraint )
    {
      more = false;
    }
    else
    {
      more = true;
    }

    return;
  }
  else
  {
    i = 0;

    for ( ; ; )
    {
      if ( x[i] < x_max[i] )
      {
        x[i] = x[i] + 1;

        constraint = 0;
        for ( j = 0; j < n; j++ )
        {
          constraint = constraint + ( double ) x[j] / ( double ) x_max[j];
        }

        if ( constraint <= 1.0 )
        {
          break;
        }
      }

      x[i] = x_min[i];

      i = i + 1;

      if ( n <= i )
      {
        more = false;
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

void vector_constrained_next4 ( int n, double alpha[], int x_min[], 
  int x_max[], int x[], double q, bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VECTOR_CONSTRAINED_NEXT4 returns the "next" constrained vector.
//
//  Discussion:
//
//    This routine is similar to VECTOR_CONSTRAINED_NEXT2 and 
//    VECTOR_CONSTRAINED_NEXT3.
//
//    We consider all vectors X of dimension N whose components
//    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
//
//    We are only interested in the subset of these vectors which
//    satisfy the following constraint:
//
//      sum ( 1 <= I <= N ) ( ALPHA(I) * X(I) ) <= Q
//
//  Example:
//
//    N = 3
//    ALPHA    4.0  3.0  5.0
//    Q       20.0
//    X_MIN:   1   1   1
//    X_MAX:   5   6   4
//
//    P = 120
//
//    #  X(1)  X(2)  X(3)      Total
//
//    1    1     1     1       12.0
//    2    2     1     1       20.0
//    3    1     2     1       15.0
//    4    2     2     1       19.0
//    5    1     3     1       18.0
//    6    1     1     2       17.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components in the vector.
//
//    Input, double ALPHA[N], the coefficient vector.
//
//    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
//    values allowed in each component.
//
//    Input/output, integer X[N].  On first call (with MORE = FALSE),
//    the input value of X is not important.  On subsequent calls, the
//    input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Input, double Q, the limit on the sum.
//
//    Input/output, bool &MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  int i;
  int j;
  double total;

  if ( !more )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = x_min[j];
    }

    total = 0.0;
    for ( j = 0; j < n; j++ )
    {
      total = total + alpha[j] * ( double ) x[j];
    }

    if ( q < total )
    {
      more = false;
    }
    else
    {
      more = true;
    }

    return;
  }
  else
  {
    i = 0;

    for ( ; ; )
    {
      if ( x[i] < x_max[i] )
      {
        x[i] = x[i] + 1;

        total = 0;
        for ( j = 0; j < n; j++ )
        {
          total = total + alpha[j] * ( double ) x[j];
        }

        if ( total <= q )
        {
          break;
        }
      }

      x[i] = x_min[i];

      i = i + 1;

      if ( n <= i )
      {
        more = false;
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

void vector_constrained_next5 ( int n, int x[], int sum_min, int sum_max, 
  bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VECTOR_CONSTRAINED_NEXT5 returns the "next" constrained vector.
//
//  Discussion:
//
//    We consider all positive integer vectors of dimension N whose 
//    components satisfy SUM_MIN <= X(1:N) <= SUM_MAX.
//
//    This routine returns, one at a time, and in right-to-left
//    lexicographic order, exactly those vectors which satisfy
//    the constraint.
//
//  Example:
//
//    N = 3
//    SUM_MIN = 5
//    SUM_MAX = 6
//
//    #  X(1)  X(2)  X(3)     SUM
//
//    1    3     1     1        5
//    2    2     2     1        5
//    3    2     1     2        5
//    4    1     3     1        5
//    5    1     2     2        5
//    6    1     1     3        5
//
//    7    4     1     1        6
//    8    3     2     1        6
//    9    3     1     2        6
//   10    2     3     1        6
//   11    2     2     2        6
//   12    2     1     3        6
//   13    1     4     1        6
//   14    1     3     2        6
//   15    1     2     3        6
//   16    1     1     4        6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of components in the vector.
//
//    Input, integer SUM_MIN, SUM_MAX, the minimum and maximum sums..
//
//    Input/output, integer X(N).  On first call (with MORE = FALSE), 
//    the input value of X is not important.  On subsequent calls, the
//    input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Input/output, logical &MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  static int base = 0;
  int i;
  int j;
//
//  Initialization.
//
  if ( !more )
  {
    if ( sum_max < n )
    {
      more = false;
      return;
    }

    if ( sum_max < sum_min )
    {
      more = false;
      return;
    }

    more = true;

    base = sum_min;
    if ( base < n )
    {
      base = n;
    }
    x[0] = base - n + 1;
    for ( i = 1; i < n; i++ )
    {
      x[i] = 1;
    }
    return;
  }
//
//  Next element.
//
  else
  {
//
//  Search from the right, seeking an index I < N for which 1 < X(I).
//
    for ( i = n-2; 0 <= i; i-- )
    {
//
//  If you find such an I, decrease X(I) by 1, and add that to X(I+1).
//
      if ( 1 < x[i] )
      {
        x[i]   = x[i]   - 1;
        x[i+1] = x[i+1] + 1;
//
//  Now grab all the "excess" 1's from the entries to the right of X(I+1).
//
        for ( j = i+2; j < n; j++ )
        {
          if ( 1 < x[j] ) 
          {
            x[i+1] = x[i+1] + x[j] - 1;
            x[j] = 1;
          }
        }
        return;
      }
    }
//
//  The current vector is (1,1,1,...BASE-N+1).
//  If BASE < SUM_MAX, then increase BASE by 1, and start the new series.
//
    if ( base < sum_max )
    {
      base = base + 1;
      x[0] = base - n + 1;
      for ( i = 1; i < n; i++ )
      {
        x[i] = 1;
      }
      return;
    }
//
//  We returned the last legal vector on the previouis call.
//  The calculation is done.
//
    more = false;
  }

  return;
}
//****************************************************************************80

void vector_constrained_next6 ( int n, double alpha[], int x_min[], 
  int x_max[], int x[], double q_min, double q_max, bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VECTOR_CONSTRAINED_NEXT6 returns the "next" constrained vector.
//
//  Discussion:
//
//    This routine is similar to VECTOR_CONSTRAINED_NEXT2,
//    VECTOR_CONSTRAINED_NEXT3, and VECTOR_CONSTRAINED_NEXT4.
//
//    We consider all vectors X of dimension N whose components
//    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
//
//    We are only interested in the subset of these vectors which
//    satisfy the following constraint:
//
//      Q_MIN <= sum ( 1 <= I <= N ) ALPHA(I) * X(I) <= Q_MAX
//
//    This routine returns, one at a time, and in right-to-left
//    lexicographic order, exactly those vectors which satisfy
//    the constraint.
//
//  Example:
//
//    N = 3
//    ALPHA    4.0  3.0  5.0
//    Q_MIN   16.0
//    Q_MAX   20.0
//    X_MIN:   1   1   1
//    X_MAX:   5   6   4
//
//    #  X(1)  X(2)  X(3)     Total
//
//    1    2     1     1       20.0
//    2    2     2     1       19.0
//    3    1     3     1       18.0
//    4    1     1     2       17.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components in the vector.
//
//    Input, double ALPHA[N], the coefficient vector.
//
//    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
//    values allowed in each component.
//
//    Input/output, int X[N].  On first call (with MORE = FALSE),
//    the input value of X is not important.  On subsequent calls, the
//    input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Input, double Q_MIN, Q_MAX, the lower and upper
//    limits on the sum.
//
//    Input/output, bool &MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  int i;
  int j;
  double total;

  if ( !more )
  {
    more = true;
    for ( i = 0; i < n; i++ )
    {
      x[i] = x_min[i];
    }

    total = 0.0;
    for ( i = 0; i < n; i++ )
    {
      total = total + alpha[i] * ( double ) ( x[i] );
    }

    if ( q_min <= total && total <= q_max )
    {
      return;
    }
  }

  for ( ; ; )
  {
    j = n - 1;

    for ( ; ; )
    {
      if ( x[j] < x_max[j] )
      {
        break;
      }

      if ( j <= 0 )
      {
        more = false;
        return;
      }
      j = j - 1;
    }

    x[j] = x[j] + 1;
    for ( i = j + 1; i < n; i++ )
    {
      x[i] = x_min[i];
    }

    total = 0.0;
    for ( i = 0; i < n; i++ )
    {
      total = total + alpha[i] * ( double ) ( x[i] );
    }

    if ( q_min <= total && total <= q_max )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void vector_constrained_next7 ( int n, double level_weight[], int x_max[], 
  int x[], double q_min, double q_max, bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VECTOR_CONSTRAINED_NEXT7 returns the "next" constrained vector.
//
//  Discussion:
//
//    We consider vectors X of dimension N satisfying:
//
//      0 <= X(1:N) <= X_MAX(1:N).
//
//    We are only interested in the subset of these vectors which
//    satisfy the following constraint:
//
//      Q_MIN < sum ( 1 <= I <= N ) LEVEL_WEIGHT(I) * X(I) <= Q_MAX
//
//    This routine returns, one at a time, and in right-to-left
//    lexicographic order, exactly those vectors which satisfy
//    the constraint.
//
//  Example:
//
//    N = 3
//    LEVEL_WEIGHT    4.0  3.0  5.0
//    Q_MIN   16.0
//    Q_MAX   20.0
//    X_MAX:   5   6   4
//
//    #  X(1)  X(2)  X(3)     Total
//
//    1    2     1     1       20.0
//    2    2     2     1       19.0
//    3    1     3     1       18.0
//    4    1     1     2       17.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components in the vector.
//
//    Input, double LEVEL_WEIGHT[N], the coefficient vector.
//
//    Input, int X_MAX[N], the maximum values allowed in each component.
//
//    Input/output, int X[N].  On first call (with MORE = FALSE),
//    the input value of X is not important.  On subsequent calls, the
//    input value of X should be the output value from the previous call.
//    On output, (with MORE = TRUE), the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors that satisfy
//    the condition.  However, on output with MORE = FALSE, the vector
//    X is meaningless, because there are no more vectors in the list.
//
//    Input, double Q_MIN, Q_MAX, the lower and upper
//    limits on the sum.
//
//    Input/output, bool &MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  int i;
  int j;
  double total;

  if ( !more )
  {
    more = true;
    for ( i = 0; i < n; i++ )
    {
      x[i] = 0;
    }

    total = 0.0;
    for ( i = 0; i < n; i++ )
    {
      total = total + level_weight[i] * ( double ) ( x[i] );
    }

    if ( q_min < total && total <= q_max )
    {
      return;
    }
  }

  for ( ; ; )
  {
    j = n - 1;

    for ( ; ; )
    {
      if ( x[j] < x_max[j] )
      {
        break;
      }

      if ( j <= 0 )
      {
        more = false;
        return;
      }
      j = j - 1;
    }

    x[j] = x[j] + 1;
    for ( i = j + 1; i < n; i++ )
    {
      x[i] = 0;
    }

    total = 0.0;
    for ( i = 0; i < n; i++ )
    {
      total = total + level_weight[i] * ( double ) ( x[i] );
    }

    if ( q_min < total && total <= q_max )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void vector_next ( int n, int x_min[], int x_max[], int x[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    VECTOR_NEXT returns the "next" integer vector between two ranges.
//
//  Discussion:
//
//    We consider all integer vectors of dimension N satisfying:
//
//      X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
//
//    This routine returns, one at a time, and in right-to-left
//    lexicographic order, all these vectors.
//
//  Example:
//
//    N = 3
//    X_MIN:   2   2   0
//    X_MAX:   4   3   1
// 
//    #  X(1)  X(2)  X(3)
//
//    1    2     2     0
//    2    3     2     0
//    3    4     2     0
//    4    2     3     0
//    5    3     3     0
//    6    4     3     0
//    7    2     2     1
//    8    3     2     1
//    9    4     2     1
//   10    2     3     1
//   11    3     3     1
//   12    4     3     1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components in the vector.
//
//    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
//    values allowed in each component.
//
//    Input/output, int X[N].  On first call, with 
//    MORE = FALSE, the input value of X is not important.  On subsequent calls,
//    the input value of X should be the output value from the previous call.
//    On output, with MORE = TRUE, the value of X will be the "next"
//    vector in the reverse lexicographical list of vectors.  However, on 
//    output with MORE = FALSE, the vector X is meaningless, because there 
//    are no more vectors in the list.
//
//    Input/output, bool &MORE.  On input, if the user has set MORE
//    FALSE, the user is requesting the initiation of a new sequence
//    of values.  If MORE is TRUE, then the user is requesting "more"
//    values in the current sequence.  On output, if MORE is TRUE,
//    then another value was found and returned in X, but if MORE is
//    FALSE, then there are no more values in the sequence, and X is
//    NOT the next value.
//
{
  int i;

  if ( !more )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = x_min[i];
    }
    more = true;
  }
  else
  {
    i = 0;

    for ( ; ; )
    {
      if ( x[i] < x_max[i] )
      {
        x[i] = x[i] + 1;
        break;
      }

      x[i] = x_min[i];

      i = i + 1;

      if ( n <= i )
      {
        more = false;
        break;
      }
    }
  }
  return;
}
//****************************************************************************80

int ytb_enum ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    YTB_ENUM enumerates the Young tableau of size N.
//
//  Discussion:
//
//    If A(N) is the number of Young tableau of size N, then A(1) = 1,
//    A(2) = 2, and
//
//    A(N) = A(N-1) + (N-1) * A(N-2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer which is to be partitioned.
//
//    Output, int YTB_ENUM, the number of Young tableau of N.
//
{
  int a1;
  int a2;
  int a3;
  int i;
  int num;

  if ( n <= 0 )
  {
    num = 0;
  }
  else if ( n == 1 )
  {
    num = 1;
  }
  else if ( n == 2 )
  {
    num = 2;
  }
  else
  {
    a2 = 1;
    a3 = 2;
    for ( i = 3; i <= n; i++ )
    {
      a1 = a2;
      a2 = a3;
      a3 = a2 + ( i - 1 ) * a1;
    }
    num = a3;
  }

  return num;
}
//****************************************************************************80

void ytb_next ( int n, int lambda[], int a[], bool &more )

//****************************************************************************80
//
//  Purpose:
//
//    YTB_NEXT computes the next Young tableau for a given shape.
//
//  Discussion:
//
//    When the routine is called with MORE = .FALSE. (the first time), and
//    if LAMBDA on this call has M parts, with M<N, then the user
//    must also make sure that LAMBDA(M+1) = 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer which is to be partitioned.
//
//    Output, int LAMBDA[N], contains a partition of N, that is,
//    the entries are positive integers that sum to N.
//
//    Output, int A[N].  A[I] is the row containing I
//    in the output tableau.
//
//    Input/output, bool &MORE.  Set MORE FALSE before first call.
//    It is reset to TRUE as the program returns a new tableau
//    on each call, until the last tableau is computed, when
//    the program also sets MORE = FALSE.
//
{
  int i;
  int ir;
  int it;
  int j;
  int k;
  int isave;

  it = n;

  if ( more )
  {
    lambda[0] = 1;
    for ( i = 1; i < n; i++ )
    {
      lambda[i] = 0;
    }

    isave = 0;

    for ( i = 2; i <= n; i++ )
    {
      lambda[a[i-1]-1] = lambda[a[i-1]-1] + 1;

      if ( a[i-1] < a[i-2] )
      {
        isave = i;
        break;
      }

    }

    if ( isave == 0 )
    {
      more = false;
      return;
    }

    it = lambda[a[isave-1]];

    for ( i = n; 1 <= i; i-- )
    {
      if ( lambda[i-1] == it )
      {
        a[isave-1] = i;
        lambda[i-1] = lambda[i-1] - 1;
        it = isave - 1;
        break;
      }

    }

  }

  k = 1;
  ir = 1;

  for ( ; ; )
  {
    if ( n < ir )
    {
      break;
    }

    if ( lambda[ir-1] != 0 )
    {
      a[k-1] = ir;
      lambda[ir-1] = lambda[ir-1] - 1;
      k = k + 1;
      ir = ir + 1;
      continue;
    }

    if ( it < k )
    {
      break;
    }

    ir = 1;

  }

  if ( n == 1 )
  {
    more = false;
    return;
  }

  for ( j = 1; j < n; j++ )
  {
    if ( a[j] < a[j-1] )
    {
      more = true;
      return;
    }
  }

  more = false;

  return;
}
//****************************************************************************80

void ytb_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    YTB_PRINT prints a Young tableau.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer that is partitioned.
//
//    Input, int A[N], describes the Young tableau.
//    A[I] is the row of the tableau on which I occurs.
//
//    Input, string TITLE, a title.
//
{
  int j;
  int *jarray;
  int row_i;
  int row_length;

  jarray = new int[n];

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  row_i = 0;

  for ( ; ; )
  {
    row_i = row_i + 1;

    row_length = 0;

    for ( j = 0; j < n; j++ )
    {
      if ( a[j] == row_i )
      {
        jarray[row_length] = j;
        row_length = row_length + 1;
      }

    }

    if ( row_length <= 0 )
    {
      break;
    }

    for ( j = 0; j < row_length; j++ )
    {
      cout << setw(6) << jarray[j]+1 << "  ";
    }
    cout << "\n";

  }

  delete [] jarray;

  return;
}
//****************************************************************************80

void ytb_random ( int n, int lambda[], int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    YTB_RANDOM selects a random Young tableau of a given shape.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 December 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer which has been partitioned.
//
//    Input, int LAMBDA[N], is a partition of N, that is,
//    N = sum ( 0 <= I < N ) LAMBDA[I].
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int A[N], the vector describing the Young tableau.
//
{
  int i;
  int ih;
  int j;
  int k;
  int m;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }

  i = 0;
  k = 0;

  for ( ; ; )
  {
    i = i + 1;

    for ( j = 0; j < lambda[i-1]; j++ )
    {
      a[j] = a[j] + 1;
      k = k + 1;
    }

    if ( n <= k )
    {
      break;
    }

  }

  for ( m = 1; m <= n; m++ )
  {

    for ( ; ; )
    {
      i = i4_uniform ( 1, a[0], seed );
      j = i4_uniform ( 1, lambda[0], seed );

      if ( i <= a[j-1] && j <= lambda[i-1] )
      {
        break;
      }
    }

    for ( ; ; )
    {
      ih = a[j-1] + lambda[i-1] - i - j;

      if ( ih == 0 )
      {
        break;
      }

      k = i4_uniform ( 1, ih, seed );

      if ( k <= lambda[i-1] - j )
      {
        j = j + k;
      }
      else
      {
        i = k - lambda[i-1] + i + j;
      }
    }

    lambda[i-1] = lambda[i-1] - 1;
    a[j-1] = a[j-1] - 1;
    a[n-m] = i;

  }

  for ( i = 0; i < n; i++ )
  {
    lambda[a[i]-1] = lambda[a[i]-1] + 1;
  }

  return;
}
