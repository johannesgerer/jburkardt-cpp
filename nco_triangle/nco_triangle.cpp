# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "nco_triangle.hpp"

//****************************************************************************80

void file_name_inc ( char *file_name )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_NAME_INC increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the input string contains no digits, a blank string is returned.
//
//    If a blank string is input, then an error condition results.
//
//  Examples:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a9to99.txt"     "a0to00.txt"  (wrap around)
//      "cat.txt"        " "           (no digits to increment)
//      " "              STOP!         (error)
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, character *FILE_NAME, (a pointer to) the character string
//    to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = s_len_trim ( file_name );

  if ( lens <= 0 )
  {
    cout << "\n";
    cout << "FILE_NAME_INC - Fatal error!\n";
    cout << "  Input file name is blank.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens-1; 0 <= i; i-- )
  {
    c = *(file_name+i);

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;
      if ( c == '9' )
      {
        c = '0';
        *(file_name+i) = c;
      }
      else
      {
        c = c + 1;
        *(file_name+i) = c;
        return;
      }
    }
  }

  if ( change == 0 )
  {
    strcpy ( file_name, " " );
  }

  return;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
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
//    I4_MIN returns the smaller of two I4's.
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
//  Formula:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//  Comments:
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Examples:
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
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
    cout << "\n";
    cout << "I4_MODP - Fatal error!\n";
    cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80*

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80*
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

int nco_triangle_degree ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_DEGREE returns the degree of an NCO rule for the triangle.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int NCO_TRIANGLE_DEGREE, the polynomial degree of exactness of
//    the rule.
//
{
  int degree;

  if ( 1 <= rule && rule <= 9 )
  {
    degree = rule - 1;
  }
  else
  {
    degree = -1;
    cout << "\n";
    cout << "NCO_TRIANGLE_DEGREE - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return degree;
}
//****************************************************************************80

int nco_triangle_order_num ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_ORDER_NUM returns the order of an NCO rule for the triangle.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int NCO_TRIANGLE_ORDER_NUM, the order (number of points)
//    of the rule.
//
{
  int order;
  int order_num;
  int *suborder;
  int suborder_num;

  suborder_num = nco_triangle_suborder_num ( rule );

  suborder = nco_triangle_suborder ( rule, suborder_num );

  order_num = 0;
  for ( order = 0; order < suborder_num; order++ )
  {
    order_num = order_num + suborder[order];
  }

  delete [] suborder;

  return order_num;
}
//****************************************************************************80

void nco_triangle_rule ( int rule, int order_num, double xy[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_RULE returns the points and weights of an NCO rule.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int ORDER_NUM, the order (number of points) of the rule.
//
//    Output, double XY[2*ORDER_NUM], the points of the rule.
//
//    Output, double W[ORDER_NUM], the weights of the rule.
//
{
  int k;
  int o;
  int s;
  int *suborder;
  int suborder_num;
  double *suborder_w;
  double *suborder_xyz;
//
//  Get the suborder information.
//
  suborder_num = nco_triangle_suborder_num ( rule );

  suborder_xyz = new double[3*suborder_num];
  suborder_w = new double[suborder_num];

  suborder = nco_triangle_suborder ( rule, suborder_num );

  nco_triangle_subrule ( rule, suborder_num, suborder_xyz, suborder_w );
//
//  Expand the suborder information to a full order rule.
//
  o = 0;

  for ( s = 0; s < suborder_num; s++ )
  {
    if ( suborder[s] == 1 )
    {
      xy[0+o*2] = suborder_xyz[0+s*3];
      xy[1+o*2] = suborder_xyz[1+s*3];
      w[o] = suborder_w[s];
      o = o + 1;
    }
    else if ( suborder[s] == 3 )
    {
      for ( k = 0; k < 3; k++ )
      {
        xy[0+o*2] = suborder_xyz [ i4_wrap(k,  0,2) + s*3 ];
        xy[1+o*2] = suborder_xyz [ i4_wrap(k+1,0,2) + s*3 ];
        w[o] = suborder_w[s];
        o = o + 1;
      }
    }
    else if ( suborder[s] == 6 )
    {
      for ( k = 0; k < 3; k++ )
      {
        xy[0+o*2] = suborder_xyz [ i4_wrap(k,  0,2) + s*3 ];
        xy[1+o*2] = suborder_xyz [ i4_wrap(k+1,0,2) + s*3 ];
        w[o] = suborder_w[s];
        o = o + 1;
      }

      for ( k = 0; k < 3; k++ )
      {
        xy[0+o*2] = suborder_xyz [ i4_wrap(k+1,0,2) + s*3 ];
        xy[1+o*2] = suborder_xyz [ i4_wrap(k,  0,2) + s*3 ];
        w[o] = suborder_w[s];
        o = o + 1;
      }
    }
    else
    {
      cout << "\n";
      cout << "NCO_TRIANGLE_RULE - Fatal error!\n;";
      cout << "  Illegal SUBORDER(" << s << ") = " << suborder[s] << "\n";
      exit ( 1 );
    }
  }

  delete [] suborder;
  delete [] suborder_xyz;
  delete [] suborder_w;

  return;
}
//****************************************************************************80

int nco_triangle_rule_num ( void )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_RULE_NUM returns the number of NCO rules available.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Output, int NCO_TRIANGLE_RULE_NUM, the number of rules available.
//
{
  int rule_num;

  rule_num = 9;

  return rule_num;
}
//****************************************************************************80

int *nco_triangle_suborder ( int rule, int suborder_num )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBORDER returns the suborders for an NCO rule.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int NCO_TRIANGLE_SUBORDER[SUBORDER_NUM],
//    the suborders of the rule.
//
{
  int *suborder;

  suborder = new int[suborder_num];

  if ( rule == 1 )
  {
    suborder[0] = 1;
  }
  else if ( rule == 2 )
  {
    suborder[0] = 3;
  }
  else if ( rule == 3 )
  {
    suborder[0] = 3;
    suborder[1] = 3;
  }
  else if ( rule == 4 )
  {
    suborder[0] = 3;
    suborder[1] = 6;
    suborder[2] = 1;
  }
  else if ( rule == 5 )
  {
    suborder[0] = 3;
    suborder[1] = 6;
    suborder[2] = 3;
    suborder[3] = 3;
  }
  else if ( rule == 6 )
  {
    suborder[0] = 3;
    suborder[1] = 6;
    suborder[2] = 6;
    suborder[3] = 3;
    suborder[4] = 3;
  }
  else if ( rule == 7 )
  {
    suborder[0] = 3;
    suborder[1] = 6;
    suborder[2] = 6;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 6;
    suborder[6] = 1;
  }
  else if ( rule == 8 )
  {
    suborder[0] = 3;
    suborder[1] = 6;
    suborder[2] = 6;
    suborder[3] = 3;
    suborder[4] = 6;
    suborder[5] = 6;
    suborder[6] = 3;
    suborder[7] = 3;
  }
  else if ( rule == 9 )
  {
    suborder[0] = 3;
    suborder[1] = 6;
    suborder[2] = 6;
    suborder[3] = 3;
    suborder[4] = 6;
    suborder[5] = 6;
    suborder[6] = 3;
    suborder[7] = 6;
    suborder[8] = 3;
    suborder[9] = 3;
  }
  else
  {
    cout << "\n";
    cout << "NCO_TRIANGLE_SUBORDER - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return suborder;
}
//****************************************************************************80

int nco_triangle_suborder_num ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBORDER_NUM returns the number of suborders for an NCO rule.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int NCO_TRIANGLE_SUBORDER_NUM, the number of suborders
//    of the rule.
//
{
  int suborder_num;

  if ( rule == 1 )
  {
    suborder_num = 1;
  }
  else if ( rule == 2 )
  {
    suborder_num = 1;
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
    suborder_num = 4;
  }
  else if ( rule == 6 )
  {
    suborder_num = 5;
  }
  else if ( rule == 7 )
  {
    suborder_num = 7;
  }
  else if ( rule == 8 )
  {
    suborder_num = 8;
  }
  else if ( rule == 9 )
  {
    suborder_num = 10;
  }
  else
  {
    suborder_num = -1;
    cout << "\n";
    cout << "NCO_TRIANGLE_SUBORDER_NUM - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return suborder_num;
}
//****************************************************************************80

void nco_triangle_subrule ( int rule, int suborder_num, double suborder_xyz[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBRULE returns a compressed NCO rule.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int i;
  int s;
  int suborder_w_d;
  int *suborder_w_n;
  int suborder_xyz_d;
  int *suborder_xyz_n;

  suborder_xyz_n = new int[3*suborder_num];
  suborder_w_n = new int[suborder_num];

  if ( rule == 1 )
  {
    nco_triangle_subrule_01 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 2 )
  {
    nco_triangle_subrule_02 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 3 )
  {
    nco_triangle_subrule_03 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 4 )
  {
    nco_triangle_subrule_04 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 5 )
  {
    nco_triangle_subrule_05 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 6 )
  {
    nco_triangle_subrule_06 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 7 )
  {
    nco_triangle_subrule_07 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 8 )
  {
    nco_triangle_subrule_08 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 9 )
  {
    nco_triangle_subrule_09 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else
  {
    cout << "\n";
    cout << "NCO_TRIANGLE_SUBRULE - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz[i+s*3] =
          ( double ) ( 1 + suborder_xyz_n[i+s*3] )
        / ( double ) ( 3 + suborder_xyz_d );
    }
  }
  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = ( double ) suborder_w_n[s] / ( double ) suborder_w_d;
  }

  delete [] suborder_w_n;
  delete [] suborder_xyz_n;

  return;
}
//****************************************************************************80

void nco_triangle_subrule_01 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBRULE_01 returns a compressed NCO rule 1.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_01[3*1] = {
      0, 0, 0
  };
  int suborder_xyz_d_01 = 0;
  int suborder_w_n_01[1] = { 1 };
  int suborder_w_d_01 = 1;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_01[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_01;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_01[s];
  }
  *suborder_w_d = suborder_w_d_01;

  return;
}
//****************************************************************************80

void nco_triangle_subrule_02 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBRULE_02 returns a compressed NCO rule 2.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_02[3*1] = {
      1, 0, 0
  };
  int suborder_xyz_d_02 = 1;
  int suborder_w_n_02[1] = { 1 };
  int suborder_w_d_02 = 3;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_02[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_02;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_02[s];
  }
  *suborder_w_d = suborder_w_d_02;

  return;
}
//****************************************************************************80

void nco_triangle_subrule_03 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBRULE_03 returns a compressed NCO rule 3.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_03[3*2] = {
      2, 0, 0,
      1, 1, 0
  };
  int suborder_xyz_d_03 = 2;
  int suborder_w_n_03[2] = { 7, -3 };
  int suborder_w_d_03 = 12;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_03[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_03;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_03[s];
  }
  *suborder_w_d = suborder_w_d_03;

  return;
}
//****************************************************************************80

void nco_triangle_subrule_04 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBRULE_04 returns a compressed NCO rule 4.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_04[3*3] = {
      3, 0, 0,
      2, 1, 0,
      1, 1, 1
  };
  int suborder_xyz_d_04 = 3;
  int suborder_w_n_04[3] = { 8, 3, -12 };
  int suborder_w_d_04 = 30;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_04[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_04;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_04[s];
  }
  *suborder_w_d = suborder_w_d_04;

  return;
}
//****************************************************************************80

void nco_triangle_subrule_05 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBRULE_05 returns a compressed NCO rule 5.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_05[3*4] = {
      4, 0, 0,
      3, 1, 0,
      2, 2, 0,
      2, 1, 1
  };
  int suborder_xyz_d_05 = 4;
  int suborder_w_n_05[4] = { 307, -316, 629, -64 };
  int suborder_w_d_05 = 720;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_05[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_05;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_05[s];
  }
  *suborder_w_d = suborder_w_d_05;

  return;
}
//****************************************************************************80

void nco_triangle_subrule_06 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBRULE_06 returns a compressed NCO rule 6.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_06[3*5] = {
      5, 0, 0,
      4, 1, 0,
      3, 2, 0,
      3, 1, 1,
      2, 2, 1
  };
  int suborder_xyz_d_06 = 5;
  int suborder_w_n_06[5] = { 71, -13, 57, -167, 113 };
  int suborder_w_d_06 = 315;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_06[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_06;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_06[s];
  }
  *suborder_w_d = suborder_w_d_06;

  return;
}
//****************************************************************************80

void nco_triangle_subrule_07 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBRULE_07 returns a compressed NCO rule 7.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_07[3*7] = {
      6, 0, 0,
      5, 1, 0,
      4, 2, 0,
      4, 1, 1,
      3, 3, 0,
      3, 2, 1,
      2, 2, 2
  };
  int suborder_xyz_d_07 = 6;
  int suborder_w_n_07[7] = { 767, -1257, 2901, 387, -3035, -915, 3509 };
  int suborder_w_d_07 = 2240;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_07[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_07;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_07[s];
  }
  *suborder_w_d = suborder_w_d_07;

  return;
}
//****************************************************************************80

void nco_triangle_subrule_08 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBRULE_08 returns a compressed NCO rule 8.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_08[3*8] = {
      7, 0, 0,
      6, 1, 0,
      5, 2, 0,
      5, 1, 1,
      4, 3, 0,
      4, 2, 1,
      3, 3, 1,
      3, 2, 2
  };
  int suborder_xyz_d_08 = 7;
  int suborder_w_n_08[8] = { 898, -662, 1573, -2522, -191, 2989, -5726, 1444 };
  int suborder_w_d_08 = 4536;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_08[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_08;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_08[s];
  }
  *suborder_w_d = suborder_w_d_08;

  return;
}
//****************************************************************************80

void nco_triangle_subrule_09 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_TRIANGLE_SUBRULE_09 returns a compressed NCO rule 9.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_09[3*10] = {
      8, 0, 0,
      7, 1, 0,
      6, 2, 0,
      6, 1, 1,
      5, 3, 0,
      5, 2, 1,
      4, 4, 0,
      4, 3, 1,
      4, 2, 2,
      3, 3, 2
  };
  int suborder_xyz_d_09 = 8;
  int suborder_w_n_09[10] = {
    1051445, -2366706, 6493915, 1818134, -9986439,-3757007, 12368047,
    478257, 10685542, -6437608 };
  int suborder_w_d_09 = 3628800;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_09[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_09;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_09[s];
  }
  *suborder_w_d = suborder_w_d_09;

  return;
}
//****************************************************************************80

double r8_huge ( void )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal double precision number,
//    and is usually defined in math.h, or sometimes in stdlib.h.
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
//    Output, double R8_HUGE, a "huge" R8 value.
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
//    R8_NINT returns the nearest integer to an R8.
//
//  Examples:
//
//        X         Value
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
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
}
//****************************************************************************80

void reference_to_physical_t3 ( double t[], int n, double ref[], double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 3 physical triangle and a point
//    (XSI,ETA) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y) in physical space.
//
//    Note that this routine may also be appropriate for an order 6
//    triangle, if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image triangle are straight and that the "midside" nodes in the
//    physical triangle are literally halfway along the sides of
//    the physical triangle.
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
//    0  1-----2
//    |
//    +--0--R--1-->
//
//  Modified:
//
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0) and
//    (0,1) respectively.
//
//    Input, int N, the number of objects to transform.
//
//    Input, double REF[2*N], points in the reference triangle.
//
//    Output, double PHY[2*N], corresponding points in the
//    physical triangle.
//
{
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = t[i+0*2] * ( 1.0 - ref[0+j*2] - ref[1+j*2] )
                 + t[i+1*2] *       + ref[0+j*2]
                 + t[i+2*2] *                    + ref[1+j*2];
    }
  }

  return;
}
//****************************************************************************80

int s_len_trim ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Modified:
//
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//****************************************************************************80

void timestamp ( void )

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
//****************************************************************************80

double triangle_area ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA computes the area of a triangle.
//
//  Discussion:
//
//    If the triangle's vertices are given in counter clockwise order,
//    the area will be positive.  If the triangle's vertices are given
//    in clockwise order, the area will be negative!
//
//    An earlier version of this routine always returned the absolute
//    value of the computed area.  I am convinced now that that is
//    a less useful result!  For instance, by returning the signed
//    area of a triangle, it is possible to easily compute the area
//    of a nonconvex polygon as the sum of the (possibly negative)
//    areas of triangles formed by node 1 and successive pairs of vertices.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA, the area of the triangle.
//
{
  double area;

  area = 0.5 * (
    t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) +
    t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) +
    t[0+2*2] * ( t[1+0*2] - t[1+1*2] ) );

  return area;
}
//****************************************************************************80

void triangle_points_plot ( char *file_name, double node_xy[], int node_show,
  int point_num, double point_xy[], int point_show )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_POINTS_PLOT plots a triangle and some points.
//
//  Modified:
//
//    04 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_NAME, the name of the output file.
//
//    Input, double NODE_XY[2*3], the coordinates of the nodes
//    of the triangle.
//
//    Input, int NODE_SHOW,
//   -1, do not show the triangle, or the nodes.
//    0, show the triangle, do not show the nodes;
//    1, show the triangle and the nodes;
//    2, show the triangle, the nodes and number them.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT_XY[2*POINT_NUM], the coordinates of the
//    points.
//
//    Input, int POINT_SHOW,
//    0, do not show the points;
//    1, show the points;
//    2, show the points and number them.
//
{
  int circle_size;
  int delta;
  int e;
  ofstream file_unit;
  int i;
  int node;
  int node_num = 3;
  int point;
  char string[40];
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;
//
//  We need to do some figuring here, so that we can determine
//  the range of the data, and hence the height and width
//  of the piece of paper.
//
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( x_max < node_xy[0+node*2] )
    {
      x_max = node_xy[0+node*2];
    }
  }
  for ( point = 0; point < point_num; point++ )
  {
    if ( x_max < point_xy[0+point*2] )
    {
      x_max = point_xy[0+point*2];
    }
  }

  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_xy[0+node*2] < x_min )
    {
      x_min = node_xy[0+node*2];
    }
  }
  for ( point = 0; point < point_num; point++ )
  {
    if ( point_xy[0+point*2] < x_min )
    {
      x_min = point_xy[0+point*2];
    }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( y_max < node_xy[1+node*2] )
    {
      y_max = node_xy[1+node*2];
    }
  }
  for ( point = 0; point < point_num; point++ )
  {
    if ( y_max < point_xy[1+point*2] )
    {
      y_max = point_xy[1+point*2];
    }
  }

  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_xy[1+node*2] < y_min )
    {
      y_min = node_xy[1+node*2];
    }
  }
  for ( point = 0; point < point_num; point++ )
  {
    if ( point_xy[1+point*2] < y_min )
    {
      y_min = point_xy[1+point*2];
    }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min )
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max      = y_ps_max - delta;
    y_ps_min      = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  file_unit.open ( file_name );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "TRIANGLE_POINTS_PLOT - Fatal error!\n";
    cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  file_unit << "%//PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: triangulation_order3_plot.C\n";
  file_unit << "%%Title: " << file_name << "\n";
  file_unit << "%%Pages: 1\n";
  file_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  file_unit << "%%Document-Fonts: Times-Roman\n";
  file_unit << "%%LanguageLevel: 1\n";
  file_unit << "%%EndComments\n";
  file_unit << "%%BeginProlog\n";
  file_unit << "/inch {72 mul} def\n";
  file_unit << "%%EndProlog\n";
  file_unit << "%%Page: 1 1\n";
  file_unit << "save\n";
  file_unit << "%\n";
  file_unit << "%  Set the RGB line color to very light gray.\n";
  file_unit << "%\n";
  file_unit << "0.900  0.900  0.900 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Draw a gray border around the page.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  file_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  file_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  file_unit << "stroke\n";
  file_unit << "%\n";
  file_unit << "%  Set the RGB color to black.\n";
  file_unit << "%\n";
  file_unit << "0.000  0.000  0.000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Set the font and its size.\n";
  file_unit << "%\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.50 inch scalefont\n";
  file_unit << "setfont\n";
  file_unit << "%\n";
  file_unit << "%  Print a title.\n";
  file_unit << "%\n";
  file_unit << "%  210  702  moveto\n";
  file_unit << "%  (Triangulation)  show\n";
  file_unit << "%\n";
  file_unit << "%  Define a clipping polygon.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  if ( 1 <= node_show )
  {
    circle_size = 5;

    file_unit << "%\n";
    file_unit << "%  Draw filled dots at the nodes.\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to blue.\n";
    file_unit << "%\n";
    file_unit << "0.000  0.150  0.750 setrgbcolor\n";
    file_unit << "%\n";

    for ( node = 0; node < 3; node++ )
    {
      x_ps = r8_nint (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (         node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = r8_nint (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  " << x_ps
                << "  " << y_ps
                << "  " << circle_size
                << "  0 360 arc closepath fill\n";
    }
  }
//
//  Label the nodes.
//
  if ( 2 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the nodes:\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to darker blue.\n";
    file_unit << "%\n";
    file_unit << "0.000  0.250  0.850 setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";
    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = r8_nint (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = r8_nint (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit  << "  " << x_ps
                 << "  " << y_ps + 5
                 << "  moveto (" << node+1 << ") show\n";
    }
  }
//
//  Draw the points.
//
  if ( point_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( point_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( point_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( point_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }

  if ( 1 <= point_show )
  {
    file_unit << "%\n";
    file_unit << "%  Draw filled dots at the points.\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to green.\n";
    file_unit << "%\n";
    file_unit << "0.150  0.750  0.000 setrgbcolor\n";
    file_unit << "%\n";

    for ( point = 0; point < point_num; point++ )
    {
      x_ps = r8_nint (
        ( ( x_max - point_xy[0+point*2]         ) * ( double ) ( x_ps_min )
        + (         point_xy[0+point*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                       - x_min ) );

      y_ps = r8_nint (
        ( ( y_max - point_xy[1+point*2]         ) * ( double ) ( y_ps_min )
        + (         point_xy[1+point*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                       - y_min ) );

      file_unit << "newpath  " << x_ps
                << "  " << y_ps
                << "  " << circle_size
                << "  0 360 arc closepath fill\n";
    }
  }
//
//  Label the points.
//
  if ( 2 <= point_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the point:\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to darker green.\n";
    file_unit << "%\n";
    file_unit << "0.250  0.850  0.000 setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";
    file_unit << "%\n";

    for ( point = 0; point < point_num; point++ )
    {
      x_ps = r8_nint (
        ( ( x_max - point_xy[0+point*2]         ) * ( double ) ( x_ps_min )
        + (       + point_xy[0+point*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                       - x_min ) );

      y_ps = r8_nint (
        ( ( y_max - point_xy[1+point*2]         ) * ( double ) ( y_ps_min )
        + (         point_xy[1+point*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                       - y_min ) );

      file_unit << "  " << x_ps
                << "  " << y_ps + 5
                << "  moveto (" << point+1 << ") show\n";
    }
  }
//
//  Draw the triangle.
//
  if ( 0 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to red.\n";
    file_unit << "%\n";
    file_unit << "0.900  0.200  0.100 setrgbcolor\n";
    file_unit << "%\n";
    file_unit << "%  Draw the triangle.\n";
    file_unit << "%\n";

    file_unit << "newpath\n";

    for ( i = 0; i <= 3; i++ )
    {
      node = i4_wrap ( i, 0, 2 );

      x_ps = ( r8_nint ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (         node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( r8_nint ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      if ( i == 0 )
      {
        file_unit << x_ps << "  " << y_ps << "  moveto\n";
      }
      else
      {
        file_unit << x_ps << "  " << y_ps << "  lineto\n";
      }
    }
    file_unit << "stroke\n";
  }

  file_unit << "%\n";
  file_unit << "restore  showpage\n";
  file_unit << "%\n";
  file_unit << "%  End of page.\n";
  file_unit << "%\n";
  file_unit << "%%Trailer\n";
  file_unit << "%%EOF\n";

  file_unit.close ( );

  return;
}
