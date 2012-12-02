# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "lyness_rule.hpp"

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
//    26 May 1999
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

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
//****************************************************************************80

int *i4vec_copy_new ( int n, int a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY_NEW copies an I4VEC to a "new" I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Output, int I4VEC_COPY_NEW[N], the copy of A1.
//
{
  int *a2;
  int i;

  a2 = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
//****************************************************************************80

int lyness_order ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    LYNESS_ORDER returns the order of a Lyness quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Lyness, Dennis Jespersen,
//    Moderate Degree Symmetric Quadrature Rules for the Triangle,
//    Journal of the Institute of Mathematics and its Applications,
//    Volume 15, Number 1, February 1975, pages 19-32.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int LYNESS_ORDER, the order of the rule.
//
{
  int order;

  if ( rule == 0 )
  {
    order = 1;
  }
  else if ( rule == 1 )
  {
    order = 3;
  }
  else if ( rule == 2 )
  {
    order = 4;
  }
  else if ( rule == 3 )
  {
    order = 4;
  }
  else if ( rule == 4 )
  {
    order = 7;
  }
  else if ( rule == 5 )
  {
    order = 6;
  }
  else if ( rule == 6 )
  {
    order = 10;
  }
  else if ( rule == 7 )
  {
    order = 9;
  }
  else if ( rule == 8 )
  {
    order = 7;
  }
  else if ( rule == 9 )
  {
    order = 10;
  }
  else if ( rule == 10 )
  {
    order = 12;
  }
  else if ( rule == 11 )
  {
    order = 16;
  }
  else if ( rule == 12 )
  {
    order = 13;
  }
  else if ( rule == 13 )
  {
    order = 13;
  }
  else if ( rule == 14 )
  {
    order = 16;
  }
  else if ( rule == 15 )
  {
    order = 16;
  }
  else if ( rule == 16 )
  {
    order = 21;
  }
  else if ( rule == 17 )
  {
    order = 16;
  }
  else if ( rule == 18 )
  {
    order = 19;
  }
  else if ( rule == 19 )
  {
    order = 22;
  }
  else if ( rule == 20 )
  {
    order = 27;
  }
  else if ( rule == 21 )
  {
    order = 28;
  }
  else
  {
    cerr << "\n";
    cerr << "LYNESS_ORDER - Fatal error!\n";
    cerr << "  Unrecognized rule index.\n";
    exit ( 1 );
  }

  return order;
}
//****************************************************************************80

int lyness_precision ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    LYNESS_PRECISION returns the precision of a Lyness quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Lyness, Dennis Jespersen,
//    Moderate Degree Symmetric Quadrature Rules for the Triangle,
//    Journal of the Institute of Mathematics and its Applications,
//    Volume 15, Number 1, February 1975, pages 19-32.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int LYNESS_PRECISION, the precision of the rule.
//
{
  int precision;

  if ( rule == 0 )
  {
    precision = 1;
  }
  else if ( rule == 1 )
  {
    precision = 2;
  }
  else if ( rule == 2 )
  {
    precision = 2;
  }
  else if ( rule == 3 )
  {
    precision = 3;
  }
  else if ( rule == 4 )
  {
    precision = 3;
  }
  else if ( rule == 5 )
  {
    precision = 4;
  }
  else if ( rule == 6 )
  {
    precision = 4;
  }
  else if ( rule == 7 )
  {
    precision = 4;
  }
  else if ( rule == 8 )
  {
    precision = 5;
  }
  else if ( rule == 9 )
  {
    precision = 5;
  }
  else if ( rule == 10 )
  {
    precision = 6;
  }
  else if ( rule == 11 )
  {
    precision = 6;
  }
  else if ( rule == 12 )
  {
    precision = 6;
  }
  else if ( rule == 13 )
  {
    precision = 7;
  }
  else if ( rule == 14 )
  {
    precision = 7;
  }
  else if ( rule == 15 )
  {
    precision = 8;
  }
  else if ( rule == 16 )
  {
    precision = 8;
  }
  else if ( rule == 17 )
  {
    precision = 8;
  }
  else if ( rule == 18 )
  {
    precision = 9;
  }
  else if ( rule == 19 )
  {
    precision = 9;
  }
  else if ( rule == 20 )
  {
    precision = 11;
  }
  else if ( rule == 21 )
  {
    precision = 11;
  }
  else
  {
    cout << "\n";
    cout << "LYNESS_PRECISION - Fatal error!\n";
    cout << "  Unrecognized rule index.\n";
    exit ( 1 );
  }

  return precision;
}
//****************************************************************************80

void lyness_rule ( int rule, int order, double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LYNESS_RULE returns the points and weights of a Lyness quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Lyness, Dennis Jespersen,
//    Moderate Degree Symmetric Quadrature Rules for the Triangle,
//    Journal of the Institute of Mathematics and its Applications,
//    Volume 15, Number 1, February 1975, pages 19-32.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int ORDER, the order of the rule.
//
//    Output, double W[ORDER], the weights.
//
//    Output, double X[2*ORDER], the points.
//
{
  int k;
  int o;
  int s;
  int *suborder;
  int suborder_num;
  double *sub_w;
  double *sub_xyz;
//
//  Get the suborder information.
//
  suborder_num = lyness_suborder_num ( rule );

  suborder = lyness_suborder ( rule, suborder_num );

  sub_xyz = new double[3*suborder_num];
  sub_w = new double[suborder_num];

  lyness_subrule ( rule, suborder_num, sub_xyz, sub_w );
//
//  Expand the suborder information to a full order rule.
//
  o = 0;

  for ( s = 0; s < suborder_num; s++ )
  {
    if ( suborder[s] == 1 )
    {
      x[0+o*2] = sub_xyz[0+s*3];
      x[1+o*2] = sub_xyz[1+s*3];
      w[o] = sub_w[s];
      o = o + 1;
    }
    else if ( suborder[s] == 3 )
    {
      for ( k = 0; k < 3; k++ )
      {
        x[0+o*2] = sub_xyz [ i4_wrap(k,  0,2) + s*3];
        x[1+o*2] = sub_xyz [ i4_wrap(k+1,0,2) + s*3];
        w[o] = sub_w[s] / 3.0;
        o = o + 1;
      }
    }
    else if ( suborder[s] == 6 )
    {
      for ( k = 0; k < 3; k++ )
      {
        x[0+o*2] = sub_xyz [ i4_wrap(k,  0,2) + s*3];
        x[1+o*2] = sub_xyz [ i4_wrap(k+1,0,2) + s*3];
        w[o] = sub_w[s] / 6.0;
        o = o + 1;
      }

      for ( k = 0; k < 3; k++ )
      {
        x[0+o*2] = sub_xyz [ i4_wrap(k+1,0,2) + s*3];
        x[1+o*2] = sub_xyz [ i4_wrap(k,  0,2) + s*3];
        w[o] = sub_w[s] / 6.0;
        o = o + 1;
      }
    }
    else
    {
      cout << "\n";
      cout << "LYNESS_RULE - Fatal error!\n";
      cout << "  Illegal SUBORDER[" << s << "] = " << suborder[2] << "\n";
      exit ( 1 );
    }
  }

  delete [] suborder;
  delete [] sub_xyz;
  delete [] sub_w;

  return;
}
//****************************************************************************80

int lyness_rule_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    LYNESS_RULE_NUM returns the number of Lyness quadrature rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Lyness, Dennis Jespersen,
//    Moderate Degree Symmetric Quadrature Rules for the Triangle,
//    Journal of the Institute of Mathematics and its Applications,
//    Volume 15, Number 1, February 1975, pages 19-32.
//
//  Parameters:
//
//    Output, int LYNESS_RULE_NUM, the number of rules.
//
{
  int rule_num;

  rule_num = 21;

  return rule_num;
}
//****************************************************************************80

int *lyness_suborder ( int rule, int suborder_num )

//****************************************************************************80
//
//  Purpose:
//
//    LYNESS_SUBORDER returns the suborders for a Lyness rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Lyness, Dennis Jespersen,
//    Moderate Degree Symmetric Quadrature Rules for the Triangle,
//    Journal of the Institute of Mathematics and its Applications,
//    Volume 15, Number 1, February 1975, pages 19-32.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int SUBORDER_NUM, the number of suborders 
//    of the rule.
//
//    Output, int LYNESS_SUBORDER[SUBORDER_NUM], the suborders 
//    of the rule.
//
{
  int *suborder;
  int suborder_00[1] = { 1 };
  int suborder_01[1] = { 3 };
  int suborder_02[2] = { 1, 3 };
  int suborder_03[2] = { 1, 3 };
  int suborder_04[3] = { 1, 3, 3 };
  int suborder_05[2] = { 3, 3 };
  int suborder_06[3] = { 1, 3, 6 };
  int suborder_07[3] = { 3, 3, 3 };
  int suborder_08[3] = { 1, 3, 3 };
  int suborder_09[4] = { 1, 3, 3, 3 };
  int suborder_10[3] = { 3, 3, 6 };
  int suborder_11[5] = { 1, 3, 3, 3, 6 };
  int suborder_12[4] = { 1, 3, 3, 6 };
  int suborder_13[4] = { 1, 3, 3, 6 };
  int suborder_14[5] = { 1, 3, 3, 3, 6 };
  int suborder_15[5] = { 1, 3, 3, 3, 6 };
  int suborder_16[6] = { 3, 3, 3, 3, 3, 6 };
  int suborder_17[5] = { 1, 3, 3, 3, 6 };
  int suborder_18[6] = { 1, 3, 3, 3, 3, 6 };
  int suborder_19[7] = { 1, 3, 3, 3, 3, 3, 6 };
  int suborder_20[7] = { 3, 3, 3, 3, 3, 6, 6 };
  int suborder_21[8] = { 1, 3, 3, 3, 3, 3, 6, 6 };

  suborder = new int[suborder_num];

  if ( rule == 0 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_00 );
  }
  else if ( rule == 1 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_01 );
  }
  else if ( rule == 2 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_02 );
  }
  else if ( rule == 3 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_03 );
  }
  else if ( rule == 4 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_04 );
  }
  else if ( rule == 5 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_05 );
  }
  else if ( rule == 6 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_06 );
  }
  else if ( rule == 7 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_07 );
  }
  else if ( rule == 8 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_08 );
  }
  else if ( rule == 9 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_09 );
  }
  else if ( rule == 10 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_10 );
  }
  else if ( rule == 11 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_11 );
  }
  else if ( rule == 12 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_12 );
  }
  else if ( rule == 13 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_13 );
  }
  else if ( rule == 14 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_14 );
  }
  else if ( rule == 15 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_15 );
  }
  else if ( rule == 16 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_16 );
  }
  else if ( rule == 17 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_17 );
  }
  else if ( rule == 18 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_18 );
  }
  else if ( rule == 19 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_19 );
  }
  else if ( rule == 20 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_20 );
  }
  else if ( rule == 21 )
  {
    suborder = i4vec_copy_new ( suborder_num, suborder_21 );
  }
  else
  {
    cout << "\n";
    cout << "LYNESS_SUBORDER - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }
  return suborder;
}
//****************************************************************************80

int lyness_suborder_num ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    LYNESS_SUBORDER_NUM returns the number of suborders for a Lyness rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Lyness, Dennis Jespersen,
//    Moderate Degree Symmetric Quadrature Rules for the Triangle,
//    Journal of the Institute of Mathematics and its Applications,
//    Volume 15, Number 1, February 1975, pages 19-32.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int LYNESS_SUBORDER_NUM, the number of suborders 
//    of the rule.
//
{
  int suborder_num;

  if ( rule == 0 )
  {
    suborder_num = 1;
  }
  else if ( rule == 1 )
  {
    suborder_num = 1;
  }
  else if ( rule == 2 )
  {
    suborder_num = 2;
  }
  else if ( rule == 3 )
  {
    suborder_num = 2;
  }
  else if ( rule == 4 )
  {
    suborder_num = 3;
  }
  else if ( rule == 5 )
  {
    suborder_num = 2;
  }
  else if ( rule == 6 )
  {
    suborder_num = 3;
  }
  else if ( rule == 7 )
  {
    suborder_num = 3;
  }
  else if ( rule == 8 )
  {
    suborder_num = 3;
  }
  else if ( rule == 9 )
  {
    suborder_num = 4;
  }
  else if ( rule == 10 )
  {
    suborder_num = 3;
  }
  else if ( rule == 11 )
  {
    suborder_num = 5;
  }
  else if ( rule == 12 )
  {
    suborder_num = 4;
  }
  else if ( rule == 13 )
  {
    suborder_num = 4;
  }
  else if ( rule == 14 )
  {
    suborder_num = 5;
  }
  else if ( rule == 15 )
  {
    suborder_num = 5;
  }
  else if ( rule == 16 )
  {
    suborder_num = 6;
  }
  else if ( rule == 17 )
  {
    suborder_num = 5;
  }
  else if ( rule == 18 )
  {
    suborder_num = 6;
  }
  else if ( rule == 19 )
  {
    suborder_num = 7;
  }
  else if ( rule == 20 )
  {
    suborder_num = 7;
  }
  else if ( rule == 21 )
  {
    suborder_num = 8;
  }
  else
  {
    cout << "\n";
    cout << "LYNESS_SUBORDER_NUM - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return suborder_num;
}


//****************************************************************************80

void lyness_subrule ( int rule, int suborder_num, double sub_xyz[], 
  double sub_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LYNESS_SUBRULE returns a compressed Lyness rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Lyness, Dennis Jespersen,
//    Moderate Degree Symmetric Quadrature Rules for the Triangle,
//    Journal of the Institute of Mathematics and its Applications,
//    Volume 15, Number 1, February 1975, pages 19-32.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int SUB_ORDER_NUM, the number of suborders 
//    of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the
//    suborder weights.
//
{
  int s;
  double sub_w_00[1] = { 
      1.0000000000E+00 };
  double sub_w_01[1] = { 
      1.0000000000E+00 };
  double sub_w_02[2] = { 
      0.7500000000E+00, 
      0.2500000000E+00 };
  double sub_w_03[2] = { 
      - 0.5625000000E+00, 
        1.5625000000E+00 };
  double sub_w_04[3] = { 
      0.45E+00, 
      0.15E+00, 
      0.40E+00 };
  double sub_w_05[2] = { 
      3.298552309659655E-01, 
      6.701447690340345E-01 };
  double sub_w_06[3] = { 
      0.45E+00, 
    - 0.05E+00, 
      0.60E+00 };
  double sub_w_07[3] = { 
      6.16204060378000851E-02, 
      0.18592649660480146E+00,      
      0.75245309735739840E+00 };
  double sub_w_08[3] = { 
      0.22500000000000000E+00, 
      0.37781754163448150E+00,      
      0.39718245836551852E+00 };
  double sub_w_09[4] = { 
      0.25312500000000000E+00, 
      0.03333333333333333E+00, 
      0.21333333333333333E+00,      
      0.50020833333333333E+00 };
  double sub_w_10[3] = { 
      3.503588271790222E-01, 
      1.525347191106164E-01,      
      4.971064537103575E-01 };
  double sub_w_11[5] = { 
      -0.57857142857142863E+00,    
      -5.95238095238095205E-02, 
       0.16190476190476191E+00, 
       1.2190476190476192E+00,  
       0.25714285714285712E+00 };
  double sub_w_12[4] = { 
      1.527089667883523E-01, 
      2.944076042366762E-01, 
      3.887052878418766E-01, 
      1.641781411330949E-01 };
  double sub_w_13[4] = { 
      - 1.495700444677495E-01, 
        5.268457722996828E-01, 
        1.600417068265167E-01, 
        4.626825653415500E-01 };
  double sub_w_14[5] = { 
      1.763126156005252E-01, 
      1.210901532768310E-02, 
      3.499561757697094E-01, 
      3.195119754425220E-01, 
      1.421102178595603E-01 };
  double sub_w_15[5] = { 
      1.443156076777862E-01, 
      2.852749028018549E-01, 
      9.737549286959440E-02, 
      3.096521116041552E-01, 
      1.633818850466092E-01 };
  double sub_w_16[6] = { 
      1.207273935292775E-02, 
     -8.491579879151455E-01, 
      1.042367468891334E+00, 
      1.947229791412260E-01, 
      4.511852767201322E-01, 
      1.488095238095238E-01 };
  double sub_w_17[5] = { 
      - 2.834183851113958E-01, 
        2.097208857979572E-01, 
        5.127273801480265E-02, 
        6.564896469913506E-01, 
        3.659351143072855E-01 };
  double sub_w_18[6] = { 
      9.713579628279610E-02, 
      9.400410068141950E-02, 
      2.334826230143263E-01, 
      2.389432167816273E-01, 
      7.673302697609430E-02, 
      2.597012362637364E-01 };
  double sub_w_19[7] = { 
      1.133624844599192E-01, 
      1.062573789846380E-03, 
      4.803411513859279E-02, 
      2.524243006337300E-01, 
      7.819254371487040E-02, 
      2.472227459993048E-01, 
      2.597012362637364E-01 };
  double sub_w_20[7] = { 
      4.097919300803106E-02, 
      1.085536215102866E-01, 
      2.781018986881812E-03, 
      1.779689321422668E-01, 
      2.314486047444677E-01, 
      3.140226717732234E-01, 
      1.242459578348437E-01 };
  double sub_w_21[8] = { 
      8.797730116222190E-02, 
      2.623293466120857E-02, 
      1.142447159818060E-01, 
      5.656634416839376E-02, 
      2.164790926342230E-01, 
      2.079874161166116E-01, 
      4.417430269980344E-02, 
      2.463378925757316E-01 };
  double sub_xyz_00[3*1] = {
      0.3333333333E+00,  0.3333333333E+00, 0.3333333334E+00 };
  double sub_xyz_01[3*1] = {
      0.0000000000E+00,  0.5000000000E+00, 0.5000000000E+00 };
  double sub_xyz_02[3*2] = {
      0.3333333333E+00,  0.3333333333E+00, 0.3333333334E+00, 
      1.0000000000E+00,  0.0000000000E+00, 0.0000000000E+00 };
  double sub_xyz_03[3*2] = {
      0.3333333333E+00,  0.3333333333E+00, 0.3333333334E+00, 
      0.6000000000E+00,  0.2000000000E+00, 0.2000000000E+00 };
  double sub_xyz_04[3*3] = {
      0.33333333333333333E+00,  0.33333333333333333E+00, 0.33333333333333333E+00, 
      1.00000000000000000E+00,  0.00000000000000000E+00, 0.00000000000000000E+00, 
      0.00000000000000000E+00,  0.50000000000000000E+00, 0.50000000000000000E+00 };
  double sub_xyz_05[3*2] = { 
      8.168475729804585E-01, 9.157621350977073E-02, 9.15762135097707569E-02, 
      1.081030181680702E-01, 4.459484909159649E-01, 0.44594849091596489E+00 };
  double sub_xyz_06[3*3] = {
      0.33333333333333333E+00,  0.33333333333333333E+00, 0.33333333333333333E+00, 
      1.00000000000000000E+00,  0.00000000000000000E+00, 0.00000000000000000E+00, 
      0.00000000000000000E+00,  0.78867513459481281E+00, 0.21132486540518719E+00 };
  double sub_xyz_07[3*3] = {
      1.00000000000000000E+00,  0.00000000000000000E+00, 0.00000000000000000E+00, 
      0.00000000000000000E+00,  0.50000000000000000E+00, 0.50000000000000000E+00, 
      0.62283903060710999E+00,  0.18858048469644506E+00, 0.18858048469644506E+00 };
  double sub_xyz_08[3*3] = {
      0.33333333333333333E+00,  0.33333333333333333E+00, 0.33333333333333333E+00, 
      0.79742698535308720E+00,  0.10128650732345633E+00, 0.10128650732345633E+00, 
      5.97158717897698088E-02,  0.47014206410511505E+00, 0.47014206410511505E+00 };
  double sub_xyz_09[3*4] = {
      0.33333333333333333E+00,  0.33333333333333333E+00, 0.33333333333333333E+00, 
      1.00000000000000000E+00,  0.00000000000000000E+00, 0.00000000000000000E+00, 
      0.00000000000000000E+00,  0.50000000000000000E+00, 0.50000000000000000E+00, 
      0.71428571428571430E-00,  0.14285714285714285E+00, 0.14285714285714285E+00 };
  double sub_xyz_10[3*3] = {
      5.014265096581342E-01,  2.492867451709329E-01, 0.24928674517093291E+00, 
      8.738219710169965E-01,  6.308901449150177E-02, 6.30890144915016854E-02, 
      6.365024991213939E-01,  5.314504984483216E-02, 0.31035245103377396E+00 };
  double sub_xyz_11[3*5] = {
      0.33333333333333333E+00,  0.33333333333333333E+00, 0.33333333333333333E+00, 
      1.00000000000000000E+00,  0.00000000000000000E+00, 0.00000000000000000E+00, 
      0.00000000000000000E+00,  0.50000000000000000E+00, 0.50000000000000000E+00, 
      0.50000000000000000E+00,  0.25000000000000000E+00, 0.25000000000000000E+00, 
      0.00000000000000000E+00,  0.90824829046386302E+00, 9.17517095361370244E-02 };
  double sub_xyz_12[3*4] = {
      3.333333333333333E-01,  3.333333333333333E-01,  0.33333333333333333E+00, 
      5.233837209269747E-02,  4.738308139536513E-01,  0.47383081395365129E+00, 
      6.557646607383649E-01,  1.721176696308175E-01,  0.17211766963081757E+00, 
      0.000000000000000E+00,  8.653073540834571E-01,  0.13469264591654295E+00 };
  double sub_xyz_13[3*4] = {
      3.333333333333333E-01,  3.333333333333333E-01,  0.33333333333333333E+00, 
      4.793080678419067E-01,  2.603459660790466E-01,  0.26034596607904670E+00,  
      8.697397941955675E-01,  6.513010290221623E-02,  6.51301029022163108E-02, 
      6.384441885698096E-01,  4.869031542531756E-02,  0.31286549600487290E+00 };
  double sub_xyz_14[3*5] = {
      3.333333333333333E-01,  3.333333333333333E-01,  0.33333333333333333E+00, 
      1.000000000000000E+00,  0.000000000000000E+00,  0.00000000000000000E+00, 
      6.901278795524791E-01,  1.549360602237604E-01,  0.15493606022376047E+00, 
      6.169850771237593E-02,  4.691507461438120E-01,  0.46915074614381203E+00, 
      0.000000000000000E+00,  8.392991722729236E-01,  0.16070082772707639E+00 };
  double sub_xyz_15[3*5] = {
      3.333333333333333E-01,  3.333333333333333E-01,  0.33333333333333333E+00, 
      8.141482341455413E-02,  4.592925882927229E-01,  0.45929258829272301E+00, 
      8.989055433659379E-01,  5.054722831703103E-02,  5.05472283170311024E-02, 
      6.588613844964797E-01,  1.705693077517601E-01,  0.17056930775176021E+00, 
      8.394777409957211E-03,  7.284923929554041E-01,  0.26311282963463867E+00 };
  double sub_xyz_16[3*6] = {
      1.000000000000000E+00,  0.000000000000000E+00,  0.00000000000000000E+00, 
      0.000000000000000E+00,  0.500000000000000E+00,  0.50000000000000000E+00, 
      8.637211648883667E-03,  4.956813941755582E-01,  0.49568139417555818E+00, 
      8.193444849714693E-01,  9.032775751426533E-02,  9.03277575142653888E-02, 
      5.316905005853895E-01,  2.341547497073052E-01,  0.23415474970730532E+00, 
      0.000000000000000E-01,  7.236067977499790E-01,  0.27639320225002095E+00 };
  double sub_xyz_17[3*5] = {
      3.333333333333333E-01,  3.333333333333333E-01,  0.33333333333333333E+00, 
      4.666912123569507E-02,  4.766654393821525E-01,  0.47666543938215239E+00, 
      9.324563118910393E-01,  3.377184405448033E-02,  3.37718440544803877E-02, 
      4.593042216691921E-01,  2.703478891654040E-01,  0.27034788916540398E+00, 
      5.146433548666149E-02,  7.458294907672514E-01,  0.20270617374608713E+00 };
  double sub_xyz_18[3*6] = {
      3.333333333333333E-01,  3.333333333333333E-01,  0.33333333333333333E+00, 
      2.063496160252593E-02,  4.896825191987370E-01,  0.48968251919873701E+00, 
      1.258208170141290E-01,  4.370895914929355E-01,  0.43708959149293541E+00, 
      6.235929287619356E-01,  1.882035356190322E-01,  0.18820353561903219E+00, 
      9.105409732110941E-01,  4.472951339445297E-02,  4.47295133944529688E-02, 
      3.683841205473626E-02,  7.411985987844980E-01,  0.22196298916076573E+00 };
  double sub_xyz_19[3*7] = {
      3.333333333333333E-01,  3.333333333333333E-01,  0.33333333333333333E+00, 
      1.000000000000000E+00,  0.000000000000000E+00,  0.00000000000000000E+00, 
      0.000000000000000E+00,  0.500000000000000E+00,  0.50000000000000000E+00, 
      1.004413236259677E-01,  4.497793381870162E-01,  0.44977933818701610E+00, 
      9.061051136018193E-01,  4.694744319909033E-02,  4.69474431990903329E-02, 
      6.162561745251021E-01,  1.918719127374489E-01,  0.19187191273744902E+00, 
      3.683841205473626E-02,  7.411985987844980E-01,  0.22196298916076573E+00 };
  double sub_xyz_20[3*7] = {
      9.352701037774565E-01,  3.236494811127173E-02,  3.23649481112718157E-02, 
      7.612981754348137E-01,  1.193509122825931E-01,  0.11935091228259319E+00, 
    - 6.922209654151433E-02,  5.346110482707572E-01,  0.53461104827075701E+00, 
      5.933801991374367E-01,  2.033099004312816E-01,  0.20330990043128172E+00, 
      2.020613940682885E-01,  3.989693029658558E-01,  0.39896930296585570E+00, 
      5.017813831049474E-02,  5.932012134282132E-01,  0.35662064826129203E+00, 
      2.102201653616613E-02,  8.074890031597923E-01,  0.17148898030404158E+00 };
  double sub_xyz_21[3*8] = {
      3.333333333333333E-01,  3.333333333333333E-01,  0.33333333333333333E+00, 
      9.480217181434233E-01,  2.598914092828833E-02,  2.59891409282883845E-02, 
      8.114249947041546E-01,  9.428750264792270E-02,  9.42875026479226691E-02, 
      1.072644996557060E-02,  4.946367750172147E-01,  0.49463677501721470E+00, 
      5.853132347709715E-01,  2.073433826145142E-01,  0.20734338261451427E+00, 
      1.221843885990187E-01,  4.389078057004907E-01,  0.43890780570049059E+00, 
      0.000000000000000E+00,  8.588702812826364E-01,  0.14112971871736357E+00, 
      4.484167758913055E-02,  6.779376548825902E-01,  0.27722066752827923E+00 };

  if ( rule == 0 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_00, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_00, sub_w );
  }
  else if ( rule == 1 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_01, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_01, sub_w );
  }
  else if ( rule == 2 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_02, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_02, sub_w );
  }
  else if ( rule == 3 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_03, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_03, sub_w );
  }
  else if ( rule == 4 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_04, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_04, sub_w );
  }
  else if ( rule == 5 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_05, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_05, sub_w );
  }
  else if ( rule == 6 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_06, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_06, sub_w );
  }
  else if ( rule == 7 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_07, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_07, sub_w );
  }
  else if ( rule == 8 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_08, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_08, sub_w );
  }
  else if ( rule == 9 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_09, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_09, sub_w );
  }
  else if ( rule == 10 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_10, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_10, sub_w );
  }
  else if ( rule == 11 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_11, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_11, sub_w );
  }
  else if ( rule == 12 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_12, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_12, sub_w );
  }
  else if ( rule == 13 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_13, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_13, sub_w );
  }
  else if ( rule == 14 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_14, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_14, sub_w );
  }
  else if ( rule == 15 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_15, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_15, sub_w );
  }
  else if ( rule == 16 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_16, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_16, sub_w );
  }
  else if ( rule == 17 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_17, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_17, sub_w );
  }
  else if ( rule == 18 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_18, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_18, sub_w );
  }
  else if ( rule == 19 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_19, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_19, sub_w );
  }
  else if ( rule == 20 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_20, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_20, sub_w );
  }
  else if ( rule == 21 )
  {
    r8mat_copy ( 3, suborder_num, sub_xyz_21, sub_xyz );
    r8vec_copy ( suborder_num, sub_w_21, sub_w );
  }
  else
  {
    cout << "\n";
    cout << "LYNESS_SUBRULE - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return;
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
    value = + x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

void r8mat_copy ( int m, int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COPY copies one R8MAT to another.
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
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double A2[M*N], the copy of A1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return;
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
//    Input, double TABLE[M*N], the table data.
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
    return;
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

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
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
