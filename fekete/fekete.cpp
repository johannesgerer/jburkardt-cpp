# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "fekete.hpp"

//****************************************************************************80

int fekete_degree ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_DEGREE returns the degree of a Fekete rule for the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int FEKETE_DEGREE, the polynomial degree of exactness of
//    the rule.
//
{
  int degree;

  if ( rule == 1 )
  {
    degree = 3;
  }
  else if ( rule == 2 )
  {
    degree = 6;
  }
  else if ( rule == 3 )
  {
    degree = 9;
  }
  else if ( rule == 4 )
  {
    degree = 12;
  }
  else if ( rule == 5 )
  {
    degree = 12;
  }
  else if ( rule == 6 )
  {
    degree = 15;
  }
  else if ( rule == 7 )
  {
    degree = 18;
  }
  else
  {
    degree = -1;
    cout << "\n";
    cout << "FEKETE_DEGREE - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return degree;
}
//****************************************************************************80

int fekete_order_num ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_ORDER_NUM returns the order of a Fekete rule for the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int FEKETE_ORDER_NUM, the order (number of points) of the rule.
//
{
  int order;
  int order_num;
  int *suborder;
  int suborder_num;

  suborder_num = fekete_suborder_num ( rule );

  suborder = fekete_suborder ( rule, suborder_num );

  order_num = 0;
  for ( order = 0; order < suborder_num; order++ )
  {
    order_num = order_num + suborder[order];
  }

  delete [] suborder;

  return order_num;
}
//****************************************************************************80

void fekete_rule ( int rule, int order_num, double xy[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_RULE returns the points and weights of a Fekete rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
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
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
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
  suborder_num = fekete_suborder_num ( rule );

  suborder_xyz = new double[3*suborder_num];
  suborder_w = new double[suborder_num];

  suborder = fekete_suborder ( rule, suborder_num );

  fekete_subrule ( rule, suborder_num, suborder_xyz, suborder_w );
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
      cout << "FEKETE_RULE - Fatal error!\n;";
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

int fekete_rule_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_RULE_NUM returns the number of Fekete rules available.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Output, int FEKETE_RULE_NUM, the number of rules available.
//
{
  int rule_num;

  rule_num = 7;

  return rule_num;
}
//****************************************************************************80

int *fekete_suborder ( int rule, int suborder_num )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_SUBORDER returns the suborders for a Fekete rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int FEKETE_SUBORDER[SUBORDER_NUM], the suborders of the rule.
//
{
  int *suborder;

  suborder = new int[suborder_num];

  if ( rule == 1 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 6;
  }
  else if ( rule == 2 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 6;
    suborder[5] = 6;
    suborder[6] = 6;
  }
  else if ( rule == 3 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 6;
    suborder[6] = 6;
    suborder[7] = 6;
    suborder[8] = 6;
    suborder[9] = 6;
    suborder[10] = 6;
    suborder[11] = 6;
  }
  else if ( rule == 4 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 3;
    suborder[6] = 3;
    suborder[7] = 6;
    suborder[8] = 6;
    suborder[9] = 6;
    suborder[10] = 6;
    suborder[11] = 6;
    suborder[12] = 6;
    suborder[13] = 6;
    suborder[14] = 6;
    suborder[15] = 6;
    suborder[16] = 6;
    suborder[17] = 6;
    suborder[18] = 6;
  }
  else if ( rule == 5 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 3;
    suborder[6] = 3;
    suborder[7] = 3;
    suborder[8] = 3;
    suborder[9] = 3;
    suborder[10] = 3;
    suborder[11] = 6;
    suborder[12] = 6;
    suborder[13] = 6;
    suborder[14] = 6;
    suborder[15] = 6;
    suborder[16] = 6;
    suborder[17] = 6;
    suborder[18] = 6;
    suborder[19] = 6;
    suborder[20] = 6;
  }
  else if ( rule == 6 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 3;
    suborder[6] = 3;
    suborder[7] = 3;
    suborder[8] = 3;
    suborder[9] = 3;
    suborder[10] = 6;
    suborder[11] = 6;
    suborder[12] = 6;
    suborder[13] = 6;
    suborder[14] = 6;
    suborder[15] = 6;
    suborder[16] = 6;
    suborder[17] = 6;
    suborder[18] = 6;
    suborder[19] = 6;
    suborder[20] = 6;
    suborder[21] = 6;
    suborder[22] = 6;
    suborder[23] = 6;
    suborder[24] = 6;
    suborder[25] = 6;
    suborder[26] = 6;
    suborder[27] = 6;
  }
  else if ( rule == 7 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 3;
    suborder[6] = 3;
    suborder[7] = 3;
    suborder[8] = 3;
    suborder[9] = 3;
    suborder[10] = 3;
    suborder[11] = 3;
    suborder[12] = 6;
    suborder[13] = 6;
    suborder[14] = 6;
    suborder[15] = 6;
    suborder[16] = 6;
    suborder[17] = 6;
    suborder[18] = 6;
    suborder[19] = 6;
    suborder[20] = 6;
    suborder[21] = 6;
    suborder[22] = 6;
    suborder[23] = 6;
    suborder[24] = 6;
    suborder[25] = 6;
    suborder[26] = 6;
    suborder[27] = 6;
    suborder[28] = 6;
    suborder[29] = 6;
    suborder[30] = 6;
    suborder[31] = 6;
    suborder[32] = 6;
    suborder[33] = 6;
    suborder[34] = 6;
    suborder[35] = 6;
    suborder[36] = 6;
    suborder[37] = 6;
  }
  else
  {
    cout << "\n";
    cout << "FEKETE_SUBORDER - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return suborder;
}
//****************************************************************************80

int fekete_suborder_num ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_SUBORDER_NUM returns the number of suborders for a Fekete rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int FEKETE_SUBORDER_NUM, the number of suborders of the rule.
//
{
  int suborder_num;

  if ( rule == 1 )
  {
    suborder_num = 3;
  }
  else if ( rule == 2 )
  {
    suborder_num = 7;
  }
  else if ( rule == 3 )
  {
    suborder_num = 12;
  }
  else if ( rule == 4 )
  {
    suborder_num = 19;
  }
  else if ( rule == 5 )
  {
    suborder_num = 21;
  }
  else if ( rule == 6 )
  {
    suborder_num = 28;
  }
  else if ( rule == 7 )
  {
    suborder_num = 38;
  }
  else
  {
    suborder_num = -1;
    cout << "\n";
    cout << "FEKETE_SUBORDER_NUM - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return suborder_num;
}
//****************************************************************************80

void fekete_subrule ( int rule, int suborder_num, double suborder_xyz[], 
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_SUBRULE returns a compressed Fekete rule.
//
//  Discussion:
//
//    The listed weights are twice what we want...since we want them
//    to sum to 1/2, reflecting the area of a unit triangle.  So we
//    simple halve the values before exiting this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
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
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
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
  int s;

  if ( rule == 1 )
  {
    fekete_subrule_1 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 2 )
  {
    fekete_subrule_2 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 3 )
  {
    fekete_subrule_3 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 4 )
  {
    fekete_subrule_4 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 5 )
  {
    fekete_subrule_5 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 6 )
  {
    fekete_subrule_6 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 7 )
  {
    fekete_subrule_7 ( suborder_num, suborder_xyz, suborder_w );
  }
  else
  {
    cout << "\n";
    cout << "FEKETE_SUBRULE - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = 0.5 * suborder_w[s];
  }

  return;
}
//****************************************************************************80

void fekete_subrule_1 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_SUBRULE_1 returns a compressed Fekete rule 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_1[3*3] = {
    0.3333333333, 0.3333333333, 0.3333333334, 
    1.0000000000, 0.0000000000, 0.0000000000,
    0.0000000000, 0.2763932023, 0.7236067977 };
  double suborder_w_rule_1[3] = {
    0.9000000000,
    0.0333333333,
    0.1666666667 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_1[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_1[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_1[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_1[s];
  }

  return;
}
//****************************************************************************80

void fekete_subrule_2 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_SUBRULE_2 returns a compressed Fekete rule 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_2[3*7] = {
      0.3333333333,  0.3333333333, 0.3333333334, 
      0.1063354684,  0.1063354684, 0.7873290632, 
      0.5000000000,  0.5000000000, 0.0000000000, 
      1.0000000000,  0.0000000000, 0.0000000000, 
      0.1171809171,  0.3162697959, 0.5665492870, 
      0.0000000000,  0.2655651402, 0.7344348598, 
      0.0000000000,  0.0848854223, 0.9151145777 };
  double suborder_w_rule_2[7] = {
      0.2178563571, 
      0.1104193374, 
      0.0358939762, 
      0.0004021278, 
      0.1771348660, 
      0.0272344079, 
      0.0192969460 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_2[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_2[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_2[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_2[s];
  }

  return;
}
//****************************************************************************80

void fekete_subrule_3 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_SUBRULE_3 returns a compressed Fekete rule 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_3[3*12] = {
      0.3333333333,  0.3333333333, 0.3333333334, 
      0.1704318201,  0.1704318201, 0.6591363598, 
      0.0600824712,  0.4699587644, 0.4699587644, 
      0.0489345696,  0.0489345696, 0.9021308608, 
      0.0000000000,  0.0000000000, 1.0000000000, 
      0.1784337588,  0.3252434900, 0.4963227512, 
      0.0588564879,  0.3010242110, 0.6401193011, 
      0.0551758079,  0.1543901944, 0.7904339977, 
      0.0000000000,  0.4173602935, 0.5826397065, 
      0.0000000000,  0.2610371960, 0.7389628040, 
      0.0000000000,  0.1306129092, 0.8693870908, 
      0.0000000000,  0.0402330070, 0.9597669930 };
  double suborder_w_rule_3[12] = {
      0.1096011288, 
      0.0767491008, 
      0.0646677819, 
      0.0276211659, 
      0.0013925011, 
      0.0933486453, 
      0.0619010169, 
      0.0437466450, 
      0.0114553907, 
      0.0093115568, 
      0.0078421987, 
      0.0022457501 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_3[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_3[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_3[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_3[s];
  }

  return;
}
//****************************************************************************80

void fekete_subrule_4 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_SUBRULE_4 returns a compressed Fekete rule 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_4[3*19] = {
      0.3333333333,  0.3333333333, 0.3333333334, 
      0.1988883477,  0.4005558262, 0.4005558261, 
      0.2618405201,  0.2618405201, 0.4763189598, 
      0.0807386775,  0.0807386775, 0.8385226450, 
      0.0336975736,  0.0336975736, 0.9326048528, 
      0.0000000000,  0.5000000000, 0.5000000000, 
      0.0000000000,  0.0000000000, 1.0000000000, 
      0.1089969290,  0.3837518758, 0.5072511952, 
      0.1590834479,  0.2454317980, 0.5954847541, 
      0.0887037176,  0.1697134458, 0.7415828366, 
      0.0302317829,  0.4071849276, 0.5625832895, 
      0.0748751152,  0.2874821712, 0.6376427136, 
      0.0250122615,  0.2489279690, 0.7260597695, 
      0.0262645218,  0.1206826354, 0.8530528428, 
      0.0000000000,  0.3753565349, 0.6246434651, 
      0.0000000000,  0.2585450895, 0.7414549105, 
      0.0000000000,  0.1569057655, 0.8430942345, 
      0.0000000000,  0.0768262177, 0.9231737823, 
      0.0000000000,  0.0233450767, 0.9766549233 }; 
  double suborder_w_rule_4[19] = {
      0.0626245179, 
      0.0571359417, 
      0.0545982307, 
      0.0172630326, 
      0.0142519606, 
      0.0030868485, 
      0.0004270742, 
      0.0455876390, 
      0.0496701966, 
      0.0387998322, 
      0.0335323983, 
      0.0268431561, 
      0.0237377452, 
      0.0177255972, 
      0.0043097313, 
      0.0028258057, 
      0.0030994935, 
      0.0023829062, 
      0.0009998683 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_4[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_4[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_4[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_4[s];
  }

  return;
}
//****************************************************************************80

void fekete_subrule_5 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_SUBRULE_5 returns a compressed Fekete rule 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_5[3*21] = {
      0.3333333333,  0.3333333333, 0.3333333334, 
      0.2201371125,  0.3169406831, 0.4629222044, 
      0.2201371125,  0.4629222044, 0.3169406831, 
      0.1877171129,  0.1877171129, 0.6245657742, 
      0.1403402144,  0.4298298928, 0.4298298928, 
      0.0833252778,  0.0833252778, 0.8333494444, 
      0.0664674598,  0.0252297247, 0.9083028155, 
      0.0218884020,  0.4890557990, 0.4890557990, 
      0.0252297247,  0.0664674598, 0.9083028155, 
      0.0000000000,  0.5000000000, 0.5000000000, 
      0.0000000000,  0.0000000000, 1.0000000000, 
      0.1157463404,  0.2842319093, 0.6000217503, 
      0.0672850606,  0.3971764400, 0.5355384994, 
      0.0909839531,  0.1779000668, 0.7311159801, 
      0.0318311633,  0.3025963402, 0.6655724965, 
      0.0273518579,  0.1733665506, 0.7992815915, 
      0.0000000000,  0.3753565349, 0.6246434651, 
      0.0000000000,  0.2585450895, 0.7414549105, 
      0.0000000000,  0.1569057655, 0.8430942345, 
      0.0000000000,  0.0768262177, 0.9231737823, 
      0.0000000000,  0.0233450767, 0.9766549233 };
  double suborder_w_rule_5[21] = {
      0.0485965670, 
      0.0602711576, 
      0.0602711576, 
      0.0476929767, 
      0.0453940802, 
      0.0258019417, 
      0.0122004614, 
      0.0230003812, 
      0.0122004614, 
      0.0018106475, 
     -0.0006601747, 
      0.0455413513, 
      0.0334182802, 
      0.0324896773, 
      0.0299402736, 
      0.0233477738, 
      0.0065962854, 
      0.0021485117, 
      0.0034785755, 
      0.0013990566, 
      0.0028825748 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_5[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_5[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_5[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_5[s];
  }

  return;
}
//****************************************************************************80

void fekete_subrule_6 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_SUBRULE_6 returns a compressed Fekete rule 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_6[3*28] = {
      0.3333333333,  0.3333333333, 0.3333333334, 
      0.2379370518,  0.3270403780, 0.4350225702, 
      0.3270403780,  0.2379370518, 0.4350225702, 
      0.1586078048,  0.4206960976, 0.4206960976, 
      0.2260541354,  0.2260541354, 0.5478917292, 
      0.1186657611,  0.1186657611, 0.7626684778, 
      0.0477095725,  0.4761452137, 0.4761452138, 
      0.0531173538,  0.0531173538, 0.8937652924, 
      0.0219495841,  0.0219495841, 0.9561008318, 
      0.0000000000,  0.0000000000, 1.0000000000, 
      0.1585345951,  0.3013819154, 0.5400834895, 
      0.0972525649,  0.3853507643, 0.5173966708, 
      0.0875150140,  0.2749910734, 0.6374939126, 
      0.1339547708,  0.1975591066, 0.6684861226, 
      0.0475622627,  0.3524012205, 0.6000365168, 
      0.0596194677,  0.1978887556, 0.7424917767, 
      0.0534939782,  0.1162464503, 0.8302595715, 
      0.0157189888,  0.4176001732, 0.5666808380, 
      0.0196887324,  0.2844332752, 0.6958779924, 
      0.0180698489,  0.1759511193, 0.8059790318, 
      0.0171941515,  0.0816639421, 0.9011419064, 
      0.0000000000,  0.4493368632, 0.5506631368, 
      0.0000000000,  0.3500847655, 0.6499152345, 
      0.0000000000,  0.2569702891, 0.7430297109, 
      0.0000000000,  0.1738056486, 0.8261943514, 
      0.0000000000,  0.1039958541, 0.8960041459, 
      0.0000000000,  0.0503997335, 0.9496002665,  
      0.0000000000,  0.0152159769, 0.9847840231 };
  double suborder_w_rule_6[28] = {
      0.0459710878, 
      0.0346650571, 
      0.0346650571, 
      0.0384470625, 
      0.0386013566, 
      0.0224308157, 
      0.0243531004, 
      0.0094392654, 
      0.0061105652, 
      0.0001283162, 
      0.0305412307, 
      0.0262101254, 
      0.0265367617, 
      0.0269859772, 
      0.0172635676, 
      0.0188795851, 
      0.0158224870, 
      0.0127170850, 
      0.0164489660, 
      0.0120018620, 
      0.0072268907, 
      0.0023599161, 
      0.0017624674, 
      0.0018648017, 
      0.0012975716, 
      0.0018506035, 
      0.0009919379, 
      0.0004893506 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_6[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_6[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_6[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_6[s];
  }

  return;
}
//****************************************************************************80

void fekete_subrule_7 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEKETE_SUBRULE_7 returns a compressed Fekete rule 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mark Taylor, Beth Wingate, Rachel Vincent,
//    An Algorithm for Computing Fekete Points in the Triangle,
//    SIAM Journal on Numerical Analysis,
//    Volume 38, Number 5, 2000, pages 1707-1720.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM], contains
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_7[3*38] = {
      0.3333333333,  0.3333333333, 0.3333333334, 
      0.2515553103,  0.3292984162, 0.4191462735, 
      0.3292984162,  0.2515553103, 0.4191462735, 
      0.1801930996,  0.4099034502, 0.4099034502, 
      0.2438647767,  0.2438647767, 0.5122704466, 
      0.1512564554,  0.1512564554, 0.6974870892, 
      0.0810689493,  0.4594655253, 0.4594655254, 
      0.0832757649,  0.0832757649, 0.8334484702, 
      0.0369065587,  0.0369065587, 0.9261868826, 
      0.0149574850,  0.0149574850, 0.9700850300, 
      0.0000000000,  0.5000000000, 0.5000000000, 
      0.0000000000,  0.0000000000, 1.0000000000, 
      0.1821465920,  0.3095465041, 0.5083069039, 
      0.1246901255,  0.3789288931, 0.4963809814, 
      0.1179441386,  0.2868915642, 0.5951642972, 
      0.1639418454,  0.2204868669, 0.6155712877, 
      0.0742549663,  0.3532533654, 0.5724916683, 
      0.0937816771,  0.2191980979, 0.6870202250, 
      0.0890951387,  0.1446273457, 0.7662775156, 
      0.0409065243,  0.4360543636, 0.5230391121, 
      0.0488675890,  0.2795984854, 0.6715339256, 
      0.0460342127,  0.2034211147, 0.7505446726, 
      0.0420687187,  0.1359040280, 0.8220272533, 
      0.0116377940,  0.4336892286, 0.5546729774, 
      0.0299062187,  0.3585587824, 0.6115349989, 
      0.0132313129,  0.2968103667, 0.6899583204, 
      0.0136098469,  0.2050279257, 0.7813622274, 
      0.0124869684,  0.1232146223, 0.8642984093, 
      0.0365197797,  0.0805854893, 0.8828947310, 
      0.0118637765,  0.0554881302, 0.9326480933, 
      0.0000000000,  0.4154069883, 0.5845930117, 
      0.0000000000,  0.3332475761, 0.6667524239, 
      0.0000000000,  0.2558853572, 0.7441146428, 
      0.0000000000,  0.1855459314, 0.8144540686, 
      0.0000000000,  0.1242528987, 0.8757471013, 
      0.0000000000,  0.0737697111, 0.9262302889, 
      0.0000000000,  0.0355492359, 0.9644507641, 
      0.0000000000,  0.0106941169, 0.9893058831 };
  double suborder_w_rule_7[38] = {
      0.0326079297, 
      0.0255331366, 
      0.0255331366, 
      0.0288093886, 
      0.0279490452, 
      0.0174438045, 
      0.0203594338, 
      0.0113349170, 
      0.0046614185, 
      0.0030346239, 
      0.0012508731, 
      0.0000782945, 
      0.0235716330, 
      0.0206304700, 
      0.0204028340, 
      0.0215105697, 
      0.0183482070, 
      0.0174161032, 
      0.0155972434, 
      0.0119269616, 
      0.0147074804, 
      0.0116182830, 
      0.0087639138, 
      0.0098563528, 
      0.0096342355, 
      0.0086477936, 
      0.0083868302, 
      0.0062576643, 
      0.0077839825, 
      0.0031415239, 
      0.0006513246, 
      0.0021137942, 
      0.0004393452, 
      0.0013662119, 
      0.0003331251, 
      0.0011613225, 
      0.0004342867, 
      0.0002031499 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_7[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_7[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_7[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_7[s];
  }

  return;
}
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
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a9to99.txt"     "a0to00.txt"  (wrap around)
//      "cat.txt"        " "           (no digits to increment)
//      " "              STOP!         (error)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
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
//    I4_MIN returns the smaller of two I4's.
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
//  Formula:
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
//    I4_WRAP forces an integer to lie between given limits by wrapping.
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

double r8_huge ( )

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
//  Example:
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
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
//****************************************************************************80

char *timestring ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTRING returns the current YMDHMS date as a string.
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
//    Output, char *TIMESTRING, a string containing the current YMDHMS date.
//
{
# define TIME_SIZE 40

  const struct tm *tm;
  size_t len;
  time_t now;
  char *s;

  now = time ( NULL );
  tm = localtime ( &now );

  s = new char[TIME_SIZE];

  len = strftime ( s, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  return s;
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
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
  char *date_time;
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

  date_time = timestring ( );
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
  file_unit << "%%CreationDate: " << date_time << "\n";
  delete [] date_time;
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
