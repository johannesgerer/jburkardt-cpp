# include <cmath>
# include <cstdlib>
# include <cstring>
# include <iomanip>
# include <iostream>

using namespace std;

# include "geometry.hpp"

int main ( );
void angle_box_2d_test ( );
void test001 ( );
void test002 ( );
void test0023 ( );
void test0025 ( );
void test003 ( );
void test0032 ( );
void test0035 ( );
void test004 ( );
void test0045 ( );
void test005 ( );
void test006 ( );
void test007 ( );
void test0075 ( );
void test008 ( );
void test0085 ( );
void test0087 ( );
void test009 ( );

void test010 ( );
void test011 ( );
void test012 ( );
void test0125 ( );
void test0126 ( );
void test0127 ( );
void test013 ( );
void test014 ( );
void test015 ( );
void test0155 ( );
void test0156 ( );
void test016 ( );
void test0165 ( );
void test017 ( );
void test018 ( );
void test0183 ( );
void test0185 ( );
void test019 ( );

void test020 ( );
void test0201 ( );
void test02015 ( );
void test0202 ( );
void test0203 ( );
void test02035 ( );
void test0204 ( );
void test0205 ( );
void test021 ( );
void test022 ( );
void test023 ( );
void test0232 ( );
void test0234 ( );
void test0235 ( );
void test0236 ( );
void test0238 ( );
void test024 ( );
void test0243 ( );
void test0245 ( );
void test025 ( );
void test0255 ( );
void test026 ( );
void test027 ( );
void test028 ( );
void test029 ( );

void test030 ( );
void test031 ( );
void test0315 ( );
void test032 ( );
void test0321 ( );
void test0322 ( );
void test0323 ( );
void test0325 ( );
void test0327 ( );
void test033 ( );
void test0335 ( );
void test0336 ( );
void test0337 ( );
void test034 ( );
void test0345 ( );
void test0346 ( );
void test035 ( );
void test0351 ( );
void test0352 ( );
void test038 ( );
void test0385 ( );
void test03855 ( );
void test0386 ( );
void test039 ( );

void test040 ( );
void test041 ( );
void test0415 ( );
void test0416 ( );
void test0418 ( );
void test042 ( );
void test043 ( );
void test044 ( );
void test045 ( );
void test046 ( );
void test047 ( );
void test0475 ( );
void test0477 ( );
void test0478 ( );
void test048 ( );
void test0485 ( );
void test049 ( );
void test0493 ( );
void test0495 ( );

void test050 ( );
void test051 ( );
void test052 ( );
void test053 ( );
void test054 ( );
void test055 ( );
void test056 ( );
void test057 ( );
void test058 ( );
void test059 ( );

void test060 ( );
void test061 ( );
void test0615 ( );
void test0616 ( );
void test0617 ( );
void test062 ( );
void test063 ( );
void test064 ( );
void test065 ( );
void test066 ( );
void test068 ( );
void test0685 ( );

void test0755 ( );
void test0757 ( );
void test076 ( );
void test0765 ( );
void test078 ( );
void test0782 ( );
void test0784 ( );
void test0786 ( );
void test079 ( );

void test080 ( );
void test0801 ( );
void test0803 ( );
void test0805 ( );
void polygon_solid_angle_3d_test ( );
void test081 ( );
void test082 ( );
void test0825 ( );
void test083 ( );
void test084 ( );
void test0844 ( );
void test0845 ( );
void test0846 ( );
void test085 ( );

void test170 ( );
void test171 ( );
void test1715 ( );
void test171 ( );
void test1712 ( );
void test172 ( );
void test173 ( );
void test174 ( );
void test1745 ( );
void test1746 ( );
void test175 ( );
void test176 ( );
void test177 ( );
void test178 ( );
void test1787 ( );
void test036 ( );
void test0365 ( );
void test0366 ( );
void test0367 ( );
void test0368 ( );
void test037 ( );
void test1788 ( );
void test1789 ( );
void test179 ( );

void test180 ( );
void test1805 ( );
void test181 ( );
void test182 ( );
void test183 ( );
void test1835 ( );
void test1836 ( );
void test187 ( );
void test188 ( );
void test189 ( );
void test1892 ( );
void test1893 ( );
void test1895 ( );

void test190 ( );
void test191 ( );
void test192 ( );
void test193 ( );
void test194 ( );
void test195 ( );
void test1955 ( );
void test196 ( );
void test197 ( );
void test198 ( );
void test199 ( );

void test200 ( );
void test201 ( );
void test202 ( );
void test203 ( );
void test2031 ( );
void test2032 ( );
void test20321 ( );
void test20322 ( );
void test203224 ( );
void test203225 ( );
void test20323 ( );
void test203232 ( );
void test203233 ( );
void test203234 ( );
void test203235 ( );
void test20324 ( );
void test20325 ( );
void tetrahedron_solid_angles_3d_test ( );
void test2033 ( );
void test204 ( );
void test205 ( );
void test206 ( );
void test20605 ( );
void test2061 ( );
void test2062 ( );
void test209 ( );
void test20655 ( );
void test2066 ( );
void test2094 ( );
void test2101 ( );
void test21011 ( );
void test2067 ( );
void test21015 ( );
void test2068 ( );
void test2069 ( );
void test207 ( );
void test2075 ( );
void test208 ( );

void test2102 ( );
void test2070 ( );
void test20701 ( );
void test2104 ( );
void test2105 ( );
void test211 ( );
void test2103 ( );
void test2071 ( );
void test20715 ( );
void test2095 ( );
void test2072 ( );
void test2115 ( );
void test212 ( );
void test213 ( );
void test219 ( );

void test220 ( );
void test221 ( );
void test222 ( );
void test2225 ( );
void test223 ( );
void test224 ( );
void test2245 ( );
void test225 ( );
void test226 ( );
void test227 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for GEOMETRY_PRB.
//
//  Discussion:
//
//    GEOMETRY_PRB tests the GEOMETRY library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   14 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "GEOMETRY_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the GEOMETRY library.\n";

  angle_box_2d_test ( );
  test001 ( );
  test002 ( );
  test0023 ( );
  test0025 ( );
  test003 ( );
  test0032 ( );
  test0035 ( );
  test004 ( );
  test0045 ( );
  test005 ( );
  test006 ( );
  test007 ( );
  test0075 ( );
  test008 ( );
  test0085 ( );
  test0087 ( );
  test009 ( );

  test010 ( );
  test011 ( );
  test012 ( );
  test0125 ( );
  test0126 ( );
  test0127 ( );
  test013 ( );
  test014 ( );
  test015 ( );
  test0155 ( );
  test0156 ( );
  test016 ( );
  test0165 ( );
  test017 ( );
  test018 ( );
  test0183 ( );
  test0185 ( );
  test019 ( );

  test020 ( );
  test0201 ( );
  test02015 ( );
  test0202 ( );
  test0203 ( );
  test02035 ( );
  test0204 ( );
  test0205 ( );
  test021 ( );
  test022 ( );
  test023 ( );
  test0232 ( );
  test0234 ( );
  test0235 ( );
  test0236 ( );
  test0238 ( );
  test024 ( );
  test0243 ( );
  test0245 ( );
  test025 ( );
  test0255 ( );
  test026 ( );
  test027 ( );
  test028 ( );
  test029 ( );

  test030 ( );
  test031 ( );
  test0315 ( );
  test032 ( );
  test0321 ( );
  test0322 ( );
  test0323 ( );
  test0325 ( );
  test0327 ( );
  test033 ( );
  test0335 ( );
  test0336 ( );
  test0337 ( );
  test034 ( );
  test0345 ( );
  test0346 ( );
  test035 ( );
  test0351 ( );
  test0352 ( );
  test038 ( );
  test0385 ( );
  test03855 ( );
  test0386 ( );
  test039 ( );

  test040 ( );
  test041 ( );
  test0415 ( );
  test0416 ( );
  test0418 ( );
  test042 ( );
  test043 ( );
  test044 ( );
  test045 ( );
  test046 ( );
  test047 ( );
  test0475 ( );
  test0477 ( );
  test0478 ( );
  test048 ( );
  test0485 ( );
  test049 ( );
  test0493 ( );
  test0495 ( );

  test050 ( );
  test051 ( );
  test052 ( );
  test053 ( );
  test054 ( );
  test055 ( );
  test056 ( );
  test057 ( );
  test058 ( );
  test059 ( );

  test060 ( );
  test061 ( );
  test0615 ( );
  test0616 ( );
  test0617 ( );
  test062 ( );
  test063 ( );
  test064 ( );
  test065 ( );
  test066 ( );
  test068 ( );
  test0685 ( );

  test0755 ( );
  test0757 ( );
  test076 ( );
  test0765 ( );
  test078 ( );
  test0782 ( );
  test0784 ( );
  test0786 ( );
  test079 ( );

  test080 ( );
  test0803 ( );
  test0805 ( );
  polygon_solid_angle_3d_test ( );
  test081 ( );
  test082 ( );
  test0825 ( );
  test083 ( );
  test084 ( );
  test0844 ( );
  test0845 ( );
  test0846 ( );
  test085 ( );

  test170 ( );
  test171 ( );
  test1712 ( );
  test1715 ( );
  test172 ( );
  test173 ( );
  test174 ( );
  test1745 ( );
  test1746 ( );
  test175 ( );
  test176 ( );
  test177 ( );
  test178 ( );
  test1787 ( );
  test1893 ( );
  test036 ( );
  test0365 ( );
  test0366 ( );
  test0367 ( );
  test0368 ( );
  test037 ( );
  test1788 ( );
  test1789 ( );
  test179 ( );

  test180 ( );
  test1805 ( );
  test181 ( );
  test182 ( );
  test183 ( );
  test1835 ( );
  test1836 ( );
  test187 ( );
  test188 ( );
  test189 ( );
  test1892 ( );
  test1895 ( );

  test190 ( );
  test191 ( );
  test192 ( );
  test193 ( );
  test194 ( );
  test195 ( );
  test1955 ( );
  test196 ( );
  test197 ( );
  test198 ( );
  test199 ( );

  test200 ( );
  test201 ( );
  test202 ( );
  test203 ( );
  test2031 ( );
  test2032 ( );
  test20321 ( );
  test20322 ( );
  test203224 ( );
  test203225 ( );
  test20323 ( );
  test203232 ( );
  test203233 ( );
  test203234 ( );
  test203235 ( );
  test20324 ( );
  test20325 ( );
  tetrahedron_solid_angles_3d_test ( );
  test2033 ( );
  test204 ( );
  test205 ( );
  test206 ( );
  test20605 ( );
  test2061 ( );
  test2062 ( );
  test209 ( );
  test20655 ( );
  test2066 ( );
  test2094 ( );
  test2101 ( );
  test21011 ( );
  test2067 ( );
  test21015 ( );
  test2068 ( );
  test2069 ( );
  test207 ( );
  test2075 ( );
  test208 ( );

  test2102 ( );
  test2070 ( );
  test20701 ( );
  test2104 ( );
  test2105 ( );
  test211 ( );
  test2103 ( );
  test2071 ( );
  test20715 ( );
  test2095 ( );
  test2072 ( );
  test2115 ( );
  test212 ( );
  test213 ( );
  test219 ( );

  test220 ( );
  test221 ( );
  test222 ( );
  test2225 ( );
  test223 ( );
  test224 ( );
  test2245 ( );
  test225 ( );
  test226 ( );
  test227 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "GEOMETRY_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void angle_box_2d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_BOX_2D_TEST tests ANGLE_BOX_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2015
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double dist;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];
  double p4[DIM_NUM];
  double p5[DIM_NUM];

  cout << "\n";
  cout << "ANGLE_BOX_2D_TEST\n";
  cout << "  ANGLE_BOX_2D\n";
  cout << "\n";
  cout << "  Compute P4 and P5, normal to\n";
  cout << "  line through P1 and P2, and\n";
  cout << "  line through P2 and P3,\n";
  cout << "  and DIST units from P2.\n";
//
//  These points define the lines
//    y = 0
//  and
//    y = 2x-6
//
  p1[0] = 0.0;
  p1[1] = 0.0;
  p2[0] = 3.0;
  p2[1] = 0.0;
  p3[0] = 4.0;
  p3[1] = 2.0;
  dist = 1.0;

  cout << "\n";
  cout << "  DIST " << setw(14) << dist  << "\n";
  cout << "  P1:  " << setw(14) << p1[0] << "  " << setw(14) << p1[1] << "\n";
  cout << "  P2:  " << setw(14) << p2[0] << "  " << setw(14) << p2[1] << "\n";
  cout << "  P3:  " << setw(14) << p3[0] << "  " << setw(14) << p3[1] << "\n";

  angle_box_2d ( dist, p1, p2, p3, p4, p5 );

  cout << "  P4:  " << setw(14) << p4[0] << "  " << setw(14) << p4[1] << "\n";
  cout << "  P5:  " << setw(14) << p5[0] << "  " << setw(14) << p5[1] << "\n";
//
//  These points define the lines
//    y = 0
//  and
//    y = 2x-6
//
  p1[0] = 0.0;
  p1[1] = 0.0;
  p2[0] = 3.0;
  p2[1] = 0.0;
  p3[0] = 2.0;
  p3[1] = -2.0;
  dist = 1.0;

  cout << "\n";
  cout << "  DIST " << setw(14) << dist  << "\n";
  cout << "  P1:  " << setw(14) << p1[0] << "  " << setw(14) << p1[1] << "\n";
  cout << "  P2:  " << setw(14) << p2[0] << "  " << setw(14) << p2[1] << "\n";
  cout << "  P3:  " << setw(14) << p3[0] << "  " << setw(14) << p3[1] << "\n";

  angle_box_2d ( dist, p1, p2, p3, p4, p5 );

  cout << "  P4:  " << setw(14) << p4[0] << "  " << setw(14) << p4[1] << "\n";
  cout << "  P5:  " << setw(14) << p5[0] << "  " << setw(14) << p5[1] << "\n";
//
//  By setting P1 = P2, we are asking to extend the line
//    y = 2x-6
//  from P3 to P2 through to the other side.
//
  p1[0] = 3.0;
  p1[1] = 0.0;
  p2[0] = 3.0;
  p2[1] = 0.0;
  p3[0] = 2.0;
  p3[1] = -2.0;
  dist = 1.0;

  cout << "\n";
  cout << "  DIST " << setw(14) << dist  << "\n";
  cout << "  P1:  " << setw(14) << p1[0] << "  " << setw(14) << p1[1] << "\n";
  cout << "  P2:  " << setw(14) << p2[0] << "  " << setw(14) << p2[1] << "\n";
  cout << "  P3:  " << setw(14) << p3[0] << "  " << setw(14) << p3[1] << "\n";

  angle_box_2d ( dist, p1, p2, p3, p4, p5 );

  cout << "  P4:  " << setw(14) << p4[0] << "  " << setw(14) << p4[1] << "\n";
  cout << "  P5:  " << setw(14) << p5[0] << "  " << setw(14) << p5[1] << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test001 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST001 tests ANGLE_CONTAINS_RAY_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 6

  int angle;
  int angle_num = 12;
  double angle_rad;
  bool inside;
  double p[DIM_NUM];
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];
  double pi = 3.141592653589793;
  int test;

  cout << "\n";
  cout << "TEST001\n";
  cout << "  ANGLE_CONTAINS_RAY_2D sees if a ray lies within an angle.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    if ( test == 0 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = 1.0;
      p3[1] = 1.0;
    }
    else if ( test == 1 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = 0.0;
      p3[1] = 1.0;
    }
    else if ( test == 2 )
    {
      p1[0] =  1.0;
      p1[1] = -1.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = 0.0;
      p3[1] = 1.0;
    }
    else if ( test == 3 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = -1.0;
      p3[1] =  0.0;
    }
    else if ( test == 4 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] =  0.0;
      p3[1] = -1.0;
    }
    else if ( test == 5 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] =  1.0;
      p3[1] = -0.01;
    }

    r8vec_print ( DIM_NUM, p1, "  Vertex A" );
    r8vec_print ( DIM_NUM, p2, "  Vertex B" );
    r8vec_print ( DIM_NUM, p3, "  Vertex C" );

    cout << "\n";
    cout << "       X            Y       Inside?\n";
    cout << "\n";

    for ( angle = 0; angle <= angle_num; angle++ )
    {
      angle_rad = ( double ) ( angle ) * 2.0 * pi / ( double ) angle_num;

      p[0] = cos ( angle_rad );
      p[1] = sin ( angle_rad );

      inside = angle_contains_ray_2d ( p1, p2, p3, p );

      cout << "  " << setw(12) << p[0]
           << "  " << setw(12) << p[1]
           << "  " << inside << "\n";
    }

  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test002 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST002 tests ANGLE_DEG_2D and ANGLE_RAD_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  int i;
  int angle_num = 12;
  double temp1;
  double temp2;
  double temp3;
  double thetad;
  double thetar;
  double v1[DIM_NUM];
  double v2[DIM_NUM];
  double v3[DIM_NUM];

  cout << "\n";
  cout << "TEST002\n";
  cout << "  ANGLE_DEG_2D computes an angle;\n";
  cout << "  ANGLE_RAD_ND computes an angle.\n";
  cout << "\n";
  cout << "      X           Y          Theta         atan2  ANGLE_RAD_ND, ANGLE_DEG_2D\n";
  cout << "\n";

  v1[0] = 1.0;
  v1[1] = 0.0;
  v3[0] = 0.0;
  v3[1] = 0.0;

  for ( i = 0; i <= angle_num; i++ )
  {
    thetad = ( double ) ( i ) * 360.0 / ( double ) ( angle_num );
    thetar = degrees_to_radians ( thetad );

    v2[0] = cos ( thetar );
    v2[1] = sin ( thetar );

    temp1 = radians_to_degrees ( atan2 ( v2[1], v2[0] ) );

    temp2 = angle_rad_nd ( DIM_NUM, v1, v2 );

    temp3 = angle_deg_2d ( v1, v3, v2 );

    cout << "  " << setw(10) << v2[0]
         << "  " << setw(10) << v2[1]
         << "  " << setw(10) << thetad
         << "  " << setw(10) << temp1
         << "  " << setw(10) << temp2
         << "  " << setw(10) << temp3 << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0023 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0023 tests ANGLE_HALF_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double angle_deg;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];
  double *p4;
  double r;

  cout << "\n";
  cout << "TEST0023\n";
  cout << "  ANGLE_HALF_2D computes the half angle between two rays;\n";
  cout << "  The angle is defined by the points (P1,P2,P3)\n";
  cout << "  or by the rays P2-->P3, P2-->P1\n";

  p2[0] = 5.0;
  p2[1] = 3.0;

  angle_deg = 75.0;
  r = 3.0;
  p1[0] = p2[0] + r * r8_cosd ( angle_deg );
  p1[1] = p2[1] + r * r8_sind ( angle_deg );

  angle_deg = 15.0;
  r = 2.0;
  p3[0] = p2[0] + r * r8_cosd ( angle_deg );
  p3[1] = p2[1] + r * r8_sind ( angle_deg );

  r8vec_print ( DIM_NUM, p1, "  Point P1:" );
  r8vec_print ( DIM_NUM, p2, "  Point P2:" );
  r8vec_print ( DIM_NUM, p3, "  Point P3:" );

  p4 = angle_half_2d ( p1, p2, p3 );

  r8vec_print ( DIM_NUM, p4,
    "  End point of unit ray from P2, defining half angle, P4:" );

  angle_deg = 45.0;
  r = 1.0;
  p4[0] = p2[0] + r * r8_cosd ( angle_deg );
  p4[1] = p2[1] + r * r8_sind ( angle_deg );

  r8vec_print ( DIM_NUM, p4, "  Expected value of P4:" );

  delete [] p4;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0025 tests ANGLE_RAD_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 6

  double angle;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];
  int test;

  cout << "\n";
  cout << "TEST0025\n";
  cout << "  ANGLE_RAD_2D computes the angle between two rays;\n";
  cout << "\n";

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    if ( test == 1 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = 1.0;
      p3[1] = 1.0;
    }
    else if ( test == 2 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = 0.0;
      p3[1] = 1.0;
    }
    else if ( test == 3 )
    {
      p1[0] = 1.0;
      p1[1] = -1.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = 0.0;
      p3[1] = 1.0;
    }
    else if ( test == 4 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = -1.0;
      p3[1] =  0.0;
    }
    else if ( test == 5 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] =  0.0;
      p3[1] = -1.0;
    }
    else if ( test == 6 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] =  1.0;
      p3[1] = -0.01;
    }

    r8vec_print ( DIM_NUM, p1, "  Vertex A" );
    r8vec_print ( DIM_NUM, p2, "  Vertex B" );
    r8vec_print ( DIM_NUM, p3, "  Vertex C" );

    angle = angle_rad_2d ( p1, p2, p3 );

    cout << "\n";
    cout << "  Angle = " << angle << "\n";
    cout << "\n";
  }
  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test003 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST003 tests ANGLE_RAD_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 3

  int i;
  double p1[DIM_NUM];
  double p1_test[DIM_NUM*TEST_NUM] = {
    1.0, 0.0, 0.0,
    1.0, 2.0, 3.0,
    0.0, 0.0, 1.0 };
  double p2[DIM_NUM] = { 0.0, 0.0, 0.0 };
  double p3[DIM_NUM] = { 0.0, 0.0, 1.0 };
  double temp1;
  double temp2;
  int test;

  cout << "\n";
  cout << "TEST003\n";
  cout << "  ANGLE_RAD_3D computes an angle;\n";
  cout << "\n";
  cout << "         X           Y           Z   ANGLE_RAD_3D  (Degrees)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      p1[i] = p1_test[i+test*DIM_NUM];
    }
    temp1 = angle_rad_3d ( p1, p2, p3 );
    temp2 = radians_to_degrees ( temp1 );

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(10) << p1[i];
    }
    cout << "  " << setw(10) << temp1
         << "  " << setw(10) << temp2 << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0032 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0032 tests ANGLE_TURN_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 13

  double p1[DIM_NUM];
  double p2[DIM_NUM] = { 0.0, 0.0 };
  double p3[DIM_NUM] = { 1.0, 0.0 };
  double pi = 3.141592653589793;
  int test;
  double theta;
  double theta_degrees;
  double turn;

  cout << "\n";
  cout << "TEST0032\n";
  cout << "  ANGLE_TURN_2D computes the turning angle\n";
  cout << "  defined by the line segments [P1,P2] and [P2,P3].\n";

  cout << "\n";
  cout << "  Our three points are:\n";
  cout << "\n";
  cout << "    P1 = (C,S)\n";
  cout << "    P2 = (0,0)\n";
  cout << "    P3 = (1,0)\n";
  cout << "\n";
  cout << "  C = cosine ( theta ), S = sine ( theta ).\n";
  cout << "\n";
  cout << "  Test  Theta    Turn\n";
  cout << "\n";

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    theta = 2.0 * pi * ( double ) ( test - 1 )
      / ( double ) ( TEST_NUM - 1 );

    theta_degrees = 360.0 * ( double ) ( test - 1 )
      / ( double ) ( TEST_NUM - 1 );

    p1[0] = cos ( theta );
    p1[1] = sin ( theta );

    turn = angle_turn_2d ( p1, p2, p3 );

    cout << "  " << setw(4)  << test
         << "  " << setw(5)  << theta_degrees
         << "  " << setw(14) << turn << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0035 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0035 tests ANNULUS_SECTOR_CENTROID_2D.
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
{
# define DIM_NUM 2

  double *centroid;
  double pc[DIM_NUM] = { 5.0, 3.0 };
  double r1 = 2.0;
  double r2 = 3.0;
  double theta1;
  double theta2;

  theta1 = degrees_to_radians ( 30.0 );
  theta2 = degrees_to_radians ( 60.0 );

  cout << "\n";
  cout << "TEST0035\n";
  cout << "  ANNULUS_SECTOR_CENTROID_2D computes the centroid of a\n";
  cout << "    circular annulus.\n";
  cout << "\n";
  cout << "  The circle has center        " << pc[0]
       << "  " << pc[1] << "\n";
  cout << "  The inner radius is R1 =     " << r1     << "\n";
  cout << "  The outer radius is R2 =     " << r2     << "\n";
  cout << "  The first angle is THETA1 =  " << theta1 << "\n";
  cout << "  The second angle is THETA2 = " << theta2 << "\n";

  centroid = annulus_sector_centroid_2d ( pc, r1, r2, theta1, theta2 );

  cout << "\n";
  cout << "  Centroid: "
       << setw(12) << centroid[0]
       << setw(12) << centroid[1] << "\n";

  delete [] centroid;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test004 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST004 tests R8_ACOS;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 9

  double temp1;
  double temp2;
  int test;
  double x;
  double x_test[TEST_NUM] = {
    5.0, 1.2, 1.0, 0.9, 0.5, 0.0, -0.9, -1.0, -1.01 };

  cout << "\n";
  cout << "TEST004\n";
  cout << "  R8_ACOS computes an angle with a given cosine;\n";
  cout << "\n";
  cout << "  X  R8_ACOS(X) (Degrees)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];

    temp1 = r8_acos ( x );
    temp2 = radians_to_degrees ( temp1 );

    cout << "  " << setw(12) << x
         << "  " << setw(12) << temp1
         << "  " << setw(12) << temp2 << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test0045 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0045 tests R8_ASIN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 9

  double temp1;
  double temp2;
  int test;
  double x;
  double x_test[TEST_NUM] = {
    5.0, 1.2, 1.0, 0.9, 0.5, 0.0, -0.9, -1.0, -1.01 };

  cout << "\n";
  cout << "TEST0045\n";
  cout << "  R8_ASIN computes an angle with a given sine;\n";
  cout << "\n";
  cout << "  X  R8_ASIN(X)  (Degrees)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];

    temp1 = r8_asin ( x );
    temp2 = radians_to_degrees ( temp1 );

    cout << "  " << setw(12) << x
         << "  " << setw(12) << temp1
         << "  " << setw(12) << temp2 << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests R8_ATAN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 8

  double temp1;
  double temp2;
  double temp3;
  int test;
  double x;
  double x_test[TEST_NUM] = { 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0 };
  double y;
  double y_test[TEST_NUM] = { 0.0, 1.0, 2.0, 0.0, -1.0, -1.0, -1.0, -1.0 };

  cout << "\n";
  cout << "TEST005\n";
  cout << "  R8_ATAN computes an angle with a given tangent.\n";
  cout << "\n";
  cout << "  X, Y, ATAN(Y/X), ATAN2(Y,X), R8_ATAN(Y,X)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    y = y_test[test];

    if ( x != 0.0 )
    {
      temp1 = atan ( y / x );
    }
    else
    {
      temp1 = r8_huge ( );
    }

    temp2 = atan2 ( y, x );
    temp3 = r8_atan ( y, x );

    cout << "  " << setw(12) << x
         << "  " << setw(12) << y
         << "  " << setw(12) << temp1
         << "  " << setw(12) << temp2
         << "  " << setw(12) << temp3<< "\n";
  }

  cout << "\n";
  cout << "  Repeat, but display answers in degrees.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    y = y_test[test];

    if ( x != 0.0 )
    {
      temp1 = radians_to_degrees ( atan ( y / x ) );
    }
    else
    {
      temp1 = r8_huge ( );
    }

    temp2 = radians_to_degrees ( atan2 ( y, x ) );
    temp3 = radians_to_degrees ( r8_atan ( y, x ) );

    cout << "  " << setw(12) << x
         << "  " << setw(12) << y
         << "  " << setw(12) << temp1
         << "  " << setw(12) << temp2
         << "  " << setw(12) << temp3<< "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test006 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 tests BALL_UNIT_SAMPLE_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double average[DIM_NUM];
  double average_r;
  double average_theta;
  int i;
  int j;
  int sample_num = 1000;
  int seed = 123456789;
  double temp;
  double theta;
  double *x;

  cout << "\n";
  cout << "TEST006\n";
  cout << "  For the unit ball in 2 dimensions (the disk):\n";
  cout << "  BALL_UNIT_SAMPLE_2D samples;\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    x = ball_unit_sample_2d ( seed );
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout<< "  " << setw(10) << x[j];
    }
    cout << "\n";
    delete [] x;
  }

  cout << "\n";
  cout << "  Number of sample points = " << sample_num << "\n";

  r8vec_zero ( DIM_NUM, average );

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_2d ( seed );
    for ( j = 0; j < DIM_NUM; j++ )
    {
      average[j] = average[j] + x[j];
    }
    delete [] x;
  }

  for ( j = 0; j < DIM_NUM; j++ )
  {
    average[j] = average[j] / ( double ) sample_num;
  }

  cout << "\n";
  cout << "  Now average the points, which should get a value\n";
  cout << "  close to zero, and closer as N increases.\n";
  cout << "\n";
  cout << "  Average:        ";
  for ( j = 0; j < DIM_NUM; j++ )
  {
    cout << "  " << setw(10) << average[j];
  }
  cout << "\n";

  average_r = 0.0;

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_2d ( seed );
    temp = 0.0;
    for ( j = 0; j < DIM_NUM; j++ )
    {
      temp = temp + x[j] * x[j];
    }
    average_r = average_r + sqrt ( temp );
    delete [] x;
  }

  average_r = average_r / ( double ) sample_num;

  cout << "\n";
  cout << "  Now average the distance of the points from the center,\n";
  cout << "  which should be 1/sqrt(2) = "
    << 1.0 / sqrt ( 2.0 ) << "\n";
  cout << "\n";
  cout << "  Average:        " << average_r << "\n";

  average_theta = 0.0;

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_2d ( seed );
    theta = r8_atan ( x[1], x[0] );
    average_theta = average_theta + theta;
    delete [] x;
  }

  average_theta = average_theta / ( double ) sample_num;

  cout << "\n";
  cout << "  Now average the angle THETA,\n";
  cout << "  which should be PI.\n";
  cout << "\n";
  cout << "  Average:        " << average_theta << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test007 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST007 tests BALL_UNIT_SAMPLE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double average[DIM_NUM];
  double average_phi;
  double average_r;
  double average_theta;
  int i;
  int j;
  int sample_num = 1000;
  double phi;
  double r;
  int seed = 123456789;
  double theta;
  double *x;

  cout << "\n";
  cout << "TEST007\n";
  cout << "  For the unit ball in 3 dimensions:\n";
  cout << "  BALL_UNIT_SAMPLE_3D samples;\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    x = ball_unit_sample_3d ( seed );
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << "  " << setw(8) << x[j];
    }
    cout << "\n";
    delete [] x;
  }

  cout << "\n";
  cout << "  Number of sample points = " << sample_num << "\n";

  r8vec_zero ( DIM_NUM, average );

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_3d ( seed );
    for ( j = 0; j < DIM_NUM; j++ )
    {
      average[j] = average[j] + x[j];
    }
    delete [] x;
  }

  for ( j = 0; j < DIM_NUM; j++ )
  {
    average[j] = average[j] / ( double ) ( sample_num );
  }

  cout << "\n";
  cout << "  Now average the points, which should get a value\n";
  cout << "  close to zero, and closer as sample_num increases.\n";
  cout << "\n";
  cout << "  Average:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << average[i];
  }
  cout << "\n";

  seed = 123456789;

  average_r = 0.0;

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_3d ( seed );
    r = 0.0;
    for ( j = 0; j < DIM_NUM; j++ )
    {
      r = r + x[j] * x[j];
    }
    r = sqrt ( r );
    average_r = average_r + r;
    delete [] x;
  }

  average_r = average_r / ( double ) ( sample_num );

  cout << "\n";
  cout << "  Now average the distance of the points from\n";
  cout << "  the center, which should be the \n";
  cout << "  1/2^(1/n) = " <<  pow ( 0.5, 1.0 / ( double ) ( DIM_NUM ) ) << "\n";
  cout << "\n";
  cout << "  Average:        " << average_r << "\n";

  average_theta = 0.0;

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_3d ( seed );
    theta = r8_atan ( x[1], x[0] );
    average_theta = average_theta + theta;
    delete [] x;
  }

  average_theta = average_theta / ( double ) ( sample_num );

  cout << "\n";
  cout << "  Now average the angle THETA,\n";
  cout << "  which should be PI.\n";
  cout << "\n";
  cout << "  Average:        " << average_theta << "\n";

  average_phi = 0.0;

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_3d ( seed );
    r = 0.0;
    for ( j = 0; j < DIM_NUM; j++ )
    {
      r = r + x[j] * x[j];
    }
    r = sqrt ( r );
    phi = acos ( x[2] / r );
    average_phi = average_phi + phi;
    delete [] x;
  }

  average_phi = average_phi / ( double ) ( sample_num );

  cout << "\n";
  cout << "  Now average the angle PHI,\n";
  cout << "  which should be PI/2.\n";
  cout << "\n";
  cout << "  Average:        " << average_phi << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0075 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0075 tests BALL_UNIT_SAMPLE_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double average[DIM_NUM];
  double average_phi;
  double average_r;
  double average_theta;
  int i;
  int j;
  double pi = 3.141592653589793;
  int sample_num = 1000;
  double phi;
  double r;
  int seed = 123456789;
  double theta;
  double *x;

  cout << "\n";
  cout << "TEST0075\n";
  cout << "  For the unit ball in N dimensions:\n";
  cout << "  BALL_UNIT_SAMPLE_ND samples;\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    x = ball_unit_sample_nd ( DIM_NUM, seed );
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << "  " << setw(8) << x[j];
    }
    cout << "\n";
    delete [] x;
  }

  cout << "\n";
  cout << "  Number of sample points = " << sample_num << "\n";

  r8vec_zero ( DIM_NUM, average );

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_nd ( DIM_NUM, seed );
    for ( j = 0; j < DIM_NUM; j++ )
    {
      average[j] = average[j] + x[j];
    }
    delete [] x;
  }

  for ( j = 0; j < DIM_NUM; j++ )
  {
    average[j] = average[j] / ( double ) ( sample_num );
  }

  cout << "\n";
  cout << "  Now average the points, which should get a value\n";
  cout << "  close to zero, and closer as N increases.\n";
  cout << "\n";
  cout << "  Average:        ";
  for ( j = 0; j < DIM_NUM; j++ )
  {
    cout << "  " << setw(8) << average[j];
  }
  cout << "\n";

  seed = 123456789;

  average_r = 0.0;

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_nd ( DIM_NUM, seed );
    r = 0.0;
    for ( j = 0; j < DIM_NUM; j++ )
    {
      r = r + x[j] * x[j];
    }
    r = sqrt ( r );
    average_r = average_r + r;
    delete [] x;
  }

  average_r = average_r / ( double ) ( sample_num );

  cout << "\n";
  cout << "  Now average the distance of the points from\n";
  cout << "  the center, which should be the\n";
  cout << "  1/2^(1/dim_num) = " 
       << pow ( 0.5, 1.0 / ( double ) ( DIM_NUM ) ) << "\n";
  cout << "\n";
  cout << "  Average:        " << average_r << "\n";

  average_theta = 0.0;

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_nd ( DIM_NUM, seed );
    theta = r8_atan ( x[1], x[0] );
    average_theta = average_theta + theta;
    delete [] x;
  }

  average_theta = average_theta / ( double ) ( sample_num );

  cout << "\n";
  cout << "  Now average the angle THETA,\n";
  cout << "  which should be PI.\n";
  cout << "\n";
  cout << "  Average:        " << average_theta << "\n";

  average_phi = 0.0;

  for ( i = 1; i <= sample_num; i++ )
  {
    x = ball_unit_sample_nd ( DIM_NUM, seed );
    r = 0.0;
    for ( j = 0; j < DIM_NUM; j++ )
    {
      r = r + x[j] * x[j];
    }
    r = sqrt ( r );
    phi = acos ( x[2] / r );
    average_phi = average_phi + phi;
    delete [] x;
  }

  average_phi = average_phi / ( double ) ( sample_num );

  cout << "\n";
  cout << "  Now average the angle PHI,\n";
  cout << "  which should be PI/2.\n";
  cout << "\n";
  cout << "  Average:        " << average_phi << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test008 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST008 tests BASIS_MAP_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double *a = NULL;
  double c[3*3];
  int i;
  int j;
  int k;
  double u[3*3] = {
    1.0, 2.0, 3.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 2.0 };
  double v[3*3] = {
    14.0, 4.0, 4.0,
     3.0, 1.0, 0.0,
     7.0, 3.0, 2.0 };

  cout << "\n";
  cout << "TEST008\n";
  cout << "  BASIS_MAP_3D computes the linear transform A\n";
  cout << "  which maps vectors U1, U2 and U3 to vectors\n";
  cout << "  V1, V2 and V3.\n";

  r8mat_print ( 3, 3, u, "  The matrix U" );

  r8mat_print ( 3, 3, v, "  The matrix V" );

  a = basis_map_3d ( u, v );

  if ( a == NULL )
  {
    cout << "\n";
    cout << "  The matrix [ U1 | U2 | U3 ] was singular.\n";
    cout << "  No transformation was computed.\n";
    return;
  }

  r8mat_print ( 3, 3, a, "  The transformation matrix" );

  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      c[i+j*3] = 0.0;
      for ( k = 0; k < 3; k++ )
      {
        c[i+j*3] = c[i+j*3] + a[i+k*3] * u[k+j*3];
      }
    }
  }

  r8mat_print ( 3, 3, c, "  The product matrix A * [ U1 | U2 | U3 ]" );

  delete [] a;
  return;
}
//****************************************************************************80

void test0085 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0085 tests BOX_01_CONTAINS_POINT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM
# define N 46

  int i;
  int j;
  double p[2];
  double xhi =  1.2;
  double xlo = -0.3;
  double yhi =  1.4;
  double ylo = -0.1;

  cout << "\n";
  cout << "TEST0085\n";
  cout << "  BOX_01_CONTAINS_POINT_2D reports if the unit box\n";
  cout << "  contains a point.\n";
  cout << "\n";
  cout << "  We will call the function repeatedly, and draw\n";
  cout << "  a sketch of the unit square.\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    p[1] = ( ( double ) ( N - i     ) * yhi
           + ( double ) (     i - 1 ) * ylo )
           / ( double ) ( N     - 1 );

    cout << "  ";
    for ( j = 1; j <= N; j++ )
    {
      p[0] = ( ( double ) ( N - j     ) * xlo
             + ( double ) (     j - 1 ) * xhi )
             / ( double ) ( N     - 1 );

      if ( box_01_contains_point_2d ( p ) )
      {
        cout << '*';
      }
      else
      {
        cout << '-';
      }
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test0087 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0087 tests BOX_CONTAINS_POINT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 46

  int i;
  int j;
  double p[2];
  double p1[2] = { -0.1, 0.3 };
  double p2[2] = {  1.1, 0.9 };
  double xhi =  1.2;
  double xlo = -0.3;
  double yhi =  1.4;
  double ylo = -0.1;

  cout << "\n";
  cout << "TEST0087\n";
  cout << "  BOX_CONTAINS_POINT_2D reports if a box\n";
  cout << "  contains a point.\n";
  cout << "\n";
  cout << "  We will call the function repeatedly, and draw\n";
  cout << "  a sketch of the box.\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    p[1] = ( ( double ) ( N - i     ) * yhi
           + ( double ) (     i - 1 ) * ylo )
           / ( double ) ( N     - 1 );

    cout << "  ";
    for ( j = 1; j <= N; j++ )
    {
      p[0] = ( ( double ) ( N - j     ) * xlo
             + ( double ) (     j - 1 ) * xhi )
             / ( double ) ( N     - 1 );

      if ( box_contains_point_2d ( p1, p2, p ) )
      {
        cout << '*';
      }
      else
      {
        cout << '-';
      }
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test009 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST009 tests BOX_SEGMENT_CLIP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 5

  int i;
  int ival;
  double p1[DIM_NUM] = { -10.0, 10.0 };
  double p2[DIM_NUM] = {  10.0, 20.0 };
  double pa[DIM_NUM];
  double pb[DIM_NUM];
  double qa[DIM_NUM];
  double qb[DIM_NUM];
  double p1_test[DIM_NUM*TEST_NUM] = {
     1.0,  2.0,
    -3.0, 12.0,
   -20.0, 20.0,
   -20.0, 40.0,
    10.0, 40.0 };
  double p2_test[DIM_NUM*TEST_NUM] = {
     8.0, 16.0,
     5.0, 12.0,
     7.0, 20.0,
     0.0,  0.0,
    20.0, 30.0 };
  int test;

  cout << "\n";
  cout << "TEST009\n";
  cout << "  BOX_SEGMENT_CLIP_2D clips a line with respect to a box.\n";
  cout << "\n";
  cout << "  The lower left box corner is:\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p1[i];
  }
  cout << "\n";
  cout << "\n";
  cout << "  The upper right box corner is:\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p2[i];
  }
  cout << "\n";
  cout << "\n";
  cout << "  We list the points PA and PB, and then\n";
  cout << "  the clipped values.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    r8vec_copy ( DIM_NUM, p1_test+test*DIM_NUM, pa );
    r8vec_copy ( DIM_NUM, p2_test+test*DIM_NUM, pb );
    r8vec_copy ( DIM_NUM, pa, qa );
    r8vec_copy ( DIM_NUM, pb, qb );

    ival = box_segment_clip_2d ( p1, p2, qa, qb );

    cout << "\n";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << pa[i];
    }
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << pb[i];
    }
    cout << "\n";
    cout << "\n";

    if ( ival == -1 )
    {
      cout << "  Line is outside the box.\n";
    }
    else if ( ival == 0 )
    {
      cout << "  Line is inside the box.\n";
    }
    else if ( ival == 1 )
    {
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << "  " << setw(8) << qa[i];
      }
    }
    else if ( ival == 2 )
    {
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << "  " << "        ";
      }
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << "  " << setw(8) << qb[i];
      }
      cout << "\n";
    }
    else if ( ival == 3 )
      {
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << "  " << setw(8) << qa[i];
      }
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << "  " << setw(8) << qb[i];
      }
      cout << "\n";
    }

  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test010 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST010 tests BOX_RAY_INT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double p1[DIM_NUM] = { 0.0, 0.0 };
  double p2[DIM_NUM] = { 5.0, 3.0 };
  double pa[DIM_NUM];
  double pa_test[DIM_NUM*TEST_NUM] = {
    3.0, 1.0,
    4.0, 1.0,
    3.0, 1.0 };
  double pb[DIM_NUM];
  double pb_test[DIM_NUM*TEST_NUM] = {
    5.0, 5.0,
    3.0, 1.0,
    4.0, 2.0 };
  double pc_test[DIM_NUM*TEST_NUM] = {
    4.0, 3.0,
    0.0, 1.0,
    5.0, 3.0 };
  double pint[DIM_NUM];
  int test;

  cout << "\n";
  cout << "TEST010\n";
  cout << "  For a box with coordinate line sides in 2D,\n";
  cout << "  BOX_RAY_INT_2D computes the intersection of\n";
  cout << "    a shape and a ray whose origin is within\n";
  cout << "    the shape.\n";
  cout << "\n";
  cout << "  Lower left box corner:\n";
  cout << "\n";
  cout << "  " << setw(12) << p1[0]
       << "  " << setw(12) << p1[1] << "\n";
  cout << "\n";
  cout << "  Upper right box corner:\n";
  cout << "\n";
  cout << "  " << setw(12) << p2[0]
       << "  " << setw(12) << p2[1] << "\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    r8vec_copy ( DIM_NUM, pa_test+test*DIM_NUM, pa );
    r8vec_copy ( DIM_NUM, pb_test+test*DIM_NUM, pb );

    box_ray_int_2d ( p1, p2, pa, pb, pint );

    cout << "\n";
    cout << "  Origin:       " << setw(12) << pa[0]
         << "  " << setw(12) << pa[1] << "\n";
    cout << "  Point 2:      " << setw(12) << pb[0]
         << "  " << setw(12) << pb[1] << "\n";
    cout << "  Intersection: " << setw(12) << pint[0]
         << "  " << setw(12) << pint[1] << "\n";
    cout << "  Correct:      " << setw(12) << pc_test[0+test*DIM_NUM]
         << "  " << setw(12) << pc_test[1+test*DIM_NUM] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test011 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST011 tests CIRCLE_DIA2IMP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double pc[DIM_NUM];
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double r;
  double theta;

  cout << "\n";
  cout << "TEST011\n";
  cout << "  CIRCLE_DIA2IMP_2D converts a diameter to an\n";
  cout << "  implicit circle in 2D.\n";

  theta = 2.0;

  p1[0] = 2.0 + 5.0 * cos ( theta );
  p1[1] = 2.0 + 5.0 * sin ( theta );
  p2[0] = 2.0 - 5.0 * cos ( theta );
  p2[1] = 2.0 - 5.0 * sin ( theta );

  r8vec_print ( DIM_NUM, p1, "  P1:" );
  r8vec_print ( DIM_NUM, p2, "  P2:" );

  circle_dia2imp_2d ( p1, p2, &r, pc );

  circle_imp_print_2d ( r, pc, "  The implicit circle:" );

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test012 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST012 tests CIRCLE_LUNE_AREA_2D, CIRCLE_SECTOR_AREA_2D, CIRCLE_TRIANGLE_AREA_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double area1;
  double area2;
  double area3;
  double pc[DIM_NUM] = { 0.0, 0.0 };
  double pi = 3.141592653589793;
  double r = 1.0;
  int test;
  int test_num = 12;
  double theta1;
  double theta2;

  cout << "\n";
  cout << "TEST012\n";
  cout << "  CIRCLE_LUNE_AREA_2D computes the area of a\n";
  cout << "    circular lune, defined by joining the endpoints\n";
  cout << "    of a circular arc.\n";
  cout << "  CIRCLE_SECTOR_AREA_2D computes the area of a\n";
  cout << "    circular sector, defined by joining the endpoints\n";
  cout << "    of a circular arc to the center.\n";
  cout << "  CIRCLE_TRIANGLE_AREA_2D computes the signed area of a\n";
  cout << "    triangle, defined by joining the endpoints\n";
  cout << "    of a circular arc and the center.\n";
  cout << "\n";
  cout << "      R            Theta1      Theta2        Sector       Triangle     Lune\n";
  cout << "\n";

  for ( test = 0; test <= test_num; test++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( test ) * 2.0 * pi / ( double ) ( test_num );

    area1 = circle_sector_area_2d ( r, pc, theta1, theta2 );

    area2 = circle_triangle_area_2d ( r, pc, theta1, theta2 );

    area3 = circle_lune_area_2d ( r, pc, theta1, theta2 );

    cout << "  " << setw(6)  << r
         << "  " << setw(12) << theta1
         << "  " << setw(12) << theta2
         << "  " << setw(12) << area1
         << "  " << setw(12) << area2
         << "  " << setw(12) << area3 << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0125 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0125 tests CIRCLE_LUNE_AREA_2D and SPHERE_CAP_VOLUME_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double area;
  double h;
  double haver_sine;
  double pc[DIM_NUM] = { 0.0, 0.0 };
  double pi = 3.141592653589793;
  double r = 1.0;
  int test;
  int test_num = 12;
  double theta1;
  double theta2;
  double volume1;
  double volume2;

  cout << "\n";
  cout << "TEST0125\n";
  cout << "  CIRCLE_LUNE_AREA_2D computes the area of a\n";
  cout << "    circular lune, defined by joining the endpoints\n";
  cout << "    of a circular arc (THETA1,THETA2).\n";
  cout << "  SPHERE_CAP_VOLUME_2D computes the volume (area) of a\n";
  cout << "    spherical cap, defined by a plane that cuts the\n";
  cout << "    sphere to a thickness of H units.\n";
  cout << "  SPHERE_CAP_VOLUME_ND does the same operation,\n";
  cout << "    but in N dimensions.\n";
  cout << "\n";
  cout << "  The two routines should get the same results\n";
  cout << "  if THETA1, THETA2 and H correspond.\n";
  cout << "\n";
  cout << "  Using a radius R = " << r << "\n";

  cout << "\n";
  cout << "        Theta1      Theta2      H           Lune        Cap        Cap\n";
  cout << "                                            area        vol_3d     vol_nd\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    h = 2.0 * r * ( double ) ( test ) / ( double ) ( test_num );

    haver_sine = sqrt ( r * r - ( r - h ) * ( r - h ) );

    if ( h <= r )
    {
      theta2 = r8_asin ( haver_sine / r );
    }
    else
    {
      theta2 = ( pi - r8_asin ( haver_sine / r ) );
    }

    theta1 = -theta2;

    area = circle_lune_area_2d ( r, pc, theta1, theta2 );

    volume1 = sphere_cap_volume_2d ( r, h );

    volume2 = sphere_cap_volume_nd ( DIM_NUM, r, h );

    cout << "  " << setw(10) << theta1
         << "  " << setw(10) << theta2
         << "  " << setw(10) << h
         << "  " << setw(10) << area
         << "  " << setw(10) << volume1
         << "  " << setw(10) << volume2 << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0126 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0126 tests SPHERE_CAP_VOLUME_3D and SPHERE_CAP_VOLUME_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double h;
  int test;
  int test_num = 12;
  double r = 1.0;
  double volume1;
  double volume2;

  cout << "\n";
  cout << "TEST0126\n";
  cout << "  SPHERE_CAP_VOLUME_3D computes the volume of a\n";
  cout << "    spherical cap, defined by a plane that cuts the\n";
  cout << "    sphere to a thickness of H units.\n";
  cout << "  SPHERE_CAP_VOLUME_ND does the same operation,\n";
  cout << "    but in N dimensions.\n";
  cout << "\n";
  cout << "  Using a radius R = " << r << "\n";

  cout << "\n";
  cout << "        H           Cap        Cap\n";
  cout << "               volume_3d  volume_nd\n";
  cout << "\n";

  for ( test = 0; test <= test_num; test++ )
  {
    h = 2.0 * r * ( double ) ( test ) / ( double ) ( test_num );

    volume1 = sphere_cap_volume_3d ( r, h );

    volume2 = sphere_cap_volume_nd ( DIM_NUM, r, h );

    cout << "  " << setw(12) << h
         << "  " << setw(12) << volume1
         << "  " << setw(12) << volume2 << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0127 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0127 tests SPHERE_CAP_AREA_3D and SPHERE_CAP_AREA_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double area1;
  double area2;
  double h;
  int test;
  int test_num = 12;
  double r = 1.0;

  cout << "\n";
  cout << "TEST0127\n";
  cout << "  SPHERE_CAP_AREA_3D computes the volume of a\n";
  cout << "    3D spherical cap, defined by a plane that cuts the\n";
  cout << "    sphere to a thickness of H units.\n";
  cout << "  SPHERE_CAP_AREA_ND computes the volume of an\n";
  cout << "    ND spherical cap, defined by a plane that cuts the\n";
  cout << "    sphere to a thickness of H units.\n";
  cout << "\n";
  cout << "        R           H           Cap         Cap\n";
  cout << "                                area_3d     area_nd\n";
  cout << "\n";

  for ( test = 0; test <= test_num; test++ )
  {
    h = 2.0 * r * ( double ) ( test ) / ( double ) ( test_num );

    area1 = sphere_cap_area_3d ( r, h );

    area2 = sphere_cap_area_nd ( DIM_NUM, r, h );

    cout << "  " << setw(12) << r
         << "  " << setw(12) << h
         << "  " << setw(12) << area1
         << "  " << setw(12) << area2 << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test013 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST013 tests CIRCLE_LUNE_CENTROID_2D and CIRCLE_SECTOR_CENTROID_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double *centroid1;
  double *centroid2;
  double pc[DIM_NUM] = { 5.0, 3.0 };
  double pi = 3.141592653589793;
  double r = 2.0;
  int test;
  int test_num = 12;
  double theta1 = 0.0;
  double theta2;

  cout << "\n";
  cout << "TEST013\n";
  cout << "  CIRCLE_LUNE_CENTROID_2D computes the centroid of a\n";
  cout << "    circular lune, defined by joining the endpoints\n";
  cout << "    of a circular arc.\n";
  cout << "  CIRCLE_SECTOR_CENTROID_2D computes the centroid of a\n";
  cout << "    circular sector, defined by joining the endpoints\n";
  cout << "    of a circular arc to the center.\n";

  circle_imp_print_2d ( r, pc, "  The implicit circle:" );

  cout << "\n";
  cout << "  The first angle of our lune and sector is always 0.\n";
  cout << "\n";
  cout << "                         Lune                       Sector\n";
  cout << "  THETA2           X             Y             X             Y\n";
  cout << "\n";

  for ( test = 0; test <= test_num; test++ )
  {
    theta2 = ( double ) ( test ) * 2.0 * pi / ( double ) ( test_num );

    centroid1 = circle_lune_centroid_2d ( r, pc, theta1, theta2 );

    centroid2 = circle_sector_centroid_2d ( r, pc, theta1, theta2 );

    cout << "  " << setw(12) << theta2
         << "  " << setw(12) << centroid1[0]
         << "  " << setw(12) << centroid1[1]
         << "  " << setw(12) << centroid2[0]
         << "  " << setw(12) << centroid2[1] << "\n";

    delete [] centroid1;
    delete [] centroid2;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test014 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST014 tests CIRCLE_EXP_CONTAINS_POINT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM

  int inside;
  double p1[2];
  double p2[2];
  double p3[2];
  double p4[2];

  cout << "\n";
  cout << "TEST014\n";
  cout << "  CIRCLE_EXP_CONTAINS_POINT_2D determines if a\n";
  cout << "    point lies inside a circle.\n";
  cout << "\n";
  cout << "  Possible return values are:\n";
  cout << "\n";
  cout << "  -1: The point is inside the circle.\n";
  cout << "   0: The point is on the circle.\n";
  cout << "   1: The point is outside the circle\n";
  cout << "   2: Colinear data, the point is on the line.\n";
  cout << "   3: Colinear data, the point is not on the line.\n";
  cout << "   4: Two equal data points, the point is on the line.\n";
  cout << "   5: Two equal data points, the point is not on the line.\n";
  cout << "   6: All data points equal, the point is equal.\n";
  cout << "   7: All data points equal, the point is not equal.\n";
//
//  This point is inside.
//
  p1[0] = 4.0;
  p1[1] = 2.0;

  p2[0] = 1.0;
  p2[1] = 5.0;

  p3[0] = -2.0;
  p3[1] = 2.0;

  p4[0] = 2.0;
  p4[1] = 3.0;

  inside = circle_exp_contains_point_2d ( p1, p2, p3, p4 );

  cout << "\n";
  cout << "  P1 = " << p1[0] << "  " << p1[1] << "\n";
  cout << "  P2 = " << p2[0] << "  " << p2[1] << "\n";
  cout << "  P3 = " << p3[0] << "  " << p3[1] << "\n";
  cout << "  P4 = " << p4[0] << "  " << p4[1] << "\n";
  cout << "  INSIDE = " << inside << "\n";
//
//  This point is actually right on the circle.
//
  p1[0] = 4.0;
  p1[1] = 2.0;

  p2[0] = 1.0;
  p2[1] = 5.0;

  p3[0] = -2.0;
  p3[1] = 2.0;

  p4[0] = 1.0;
  p4[1] = -1.0;

  inside = circle_exp_contains_point_2d ( p1, p2, p3, p4 );

  cout << "\n";
  cout << "  P1 = " << p1[0] << "  " << p1[1] << "\n";
  cout << "  P2 = " << p2[0] << "  " << p2[1] << "\n";
  cout << "  P3 = " << p3[0] << "  " << p3[1] << "\n";
  cout << "  P4 = " << p4[0] << "  " << p4[1] << "\n";
  cout << "  INSIDE = " << inside << "\n";
//
//  This point is outside.
//
  p1[0] = 4.0;
  p1[1] = 2.0;

  p2[0] = 1.0;
  p2[1] = 5.0;

  p3[0] = -2.0;
  p3[1] = 2.0;

  p4[0] = 4.0;
  p4[1] = 6.0;

  inside = circle_exp_contains_point_2d ( p1, p2, p3, p4 );

  cout << "\n";
  cout << "  P1 = " << p1[0] << "  " << p1[1] << "\n";
  cout << "  P2 = " << p2[0] << "  " << p2[1] << "\n";
  cout << "  P3 = " << p3[0] << "  " << p3[1] << "\n";
  cout << "  P4 = " << p4[0] << "  " << p4[1] << "\n";
  cout << "  INSIDE = " << inside << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST015 tests CIRCLE_EXP2IMP_2D and TRIANGLE_CIRCUMCIRCLE_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double *p1;
  double *p2;
  double *p3;
  double pc[DIM_NUM];
  double r;
  double *t;
  int test;
  double test_t[DIM_NUM*3*TEST_NUM] = {
    4.0, 2.0,
    1.0, 5.0,
   -2.0, 2.0,
    4.0, 2.0,
    5.0, 4.0,
    6.0, 6.0,
    4.0, 2.0,
    1.0, 5.0,
    4.0, 2.0 };

  cout << "\n";
  cout << "TEST015\n";
  cout << "  CIRCLE_EXP2IMP_2D computes the radius and\n";
  cout << "    center of the circle through three points.\n";
  cout << "  TRIANGLE_CIRCUMCIRCLE_2D computes the radius and\n";
  cout << "    center of the circle through the vertices of\n";
  cout << "    a triangle.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    t = test_t+test*DIM_NUM*3;

    r8mat_transpose_print ( DIM_NUM, 3, t, "  The triangle:" );

    p1 = t;
    p2 = p1 + DIM_NUM;
    p3 = p2 + DIM_NUM;

    circle_exp2imp_2d ( p1, p2, p3, &r, pc );

    circle_imp_print_2d ( r, pc, "  The implicit circle:" );

    triangle_circumcircle_2d ( t, &r, pc );

    circle_imp_print_2d ( r, pc, "  The triangle's circumcircle:" );
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test0155 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0155 tests CIRCLE_EXP2IMP_2D.
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
{
# define DIM_NUM 2
# define TEST_NUM 13

  double curvature;
  double pc[DIM_NUM];
  double p1[DIM_NUM] = { 0.0, 0.0 };
  double p2[DIM_NUM] = { 1.0, 0.0 };
  double p3[DIM_NUM];
  double pi = 3.141592653589793;
  double r;
  double theta;
  double theta_degrees;
  int test;

  cout << "\n";
  cout << "TEST0155\n";
  cout << "  CIRCLE_EXP2IMP_2D computes the radius and\n";
  cout << "    center of the circle through three points.\n";
  cout << "\n";
  cout << "  We can use this routine to compute, for three\n";
  cout << "  points in space, the circle incident to those\n";
  cout << "  points, and hence the radius of that circle,\n";
  cout << "  and hence the \"curvature\" of those points.\n";

  cout << "\n";
  cout << "  Our three points are:\n";
  cout << "\n";
  cout << "    (0,0)\n";
  cout << "    (1,0)\n";
  cout << "    (C,S)\n";
  cout << "\n";
  cout << "  C = cosine ( theta), S = sine ( theta ).\n";
  cout << "\n";
  cout << "  Test  Theta  Curvature\n";
  cout << "\n";

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    theta = 2.0 * pi * ( double ) ( test     - 1 )
                     / ( double ) ( TEST_NUM - 1 );

    theta_degrees = 360.0 * ( double ) ( test     - 1 )
                         /  ( double ) ( TEST_NUM - 1 );

    p3[0] = cos ( theta );
    p3[1] = sin ( theta );

    circle_exp2imp_2d ( p1, p2, p3, &r, pc );

    if ( 0.0 < r )
    {
      curvature = 1.0 / r;
    }
    else
    {
      curvature = 0.0;
    }
    cout << "  " << setw(4)  << test
         << "  " << setw(5)  << theta_degrees
         << "  " << setw(14) << curvature << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0156 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0156 tests CIRCLE_EXP2IMP_2D and CIRCLE_IMP2EXP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];
  double pc1[DIM_NUM];
  double pc2[DIM_NUM];
  double r1;
  double r2;

  cout << "\n";
  cout << "TEST0156\n";
  cout << "  CIRCLE_EXP2IMP_2D converts an explicit circle\n";
  cout << "    to an implicit circle.\n";
  cout << "  CIRCLE_IMP2EXP_2D converts an implicit circle\n";
  cout << "    to an explicit circle.\n";

  pc1[0] = 10.0;
  pc1[1] =  5.0;
  r1 = 3.0;

  circle_imp_print_2d ( r1, pc1, "  The implicit circle:" );

  circle_imp2exp_2d ( r1, pc1, p1, p2, p3 );

  r8vec_print ( DIM_NUM, p1, "  P1:" );
  r8vec_print ( DIM_NUM, p2, "  P2:" );
  r8vec_print ( DIM_NUM, p3, "  P3:" );

  circle_exp2imp_2d ( p1, p2, p3, &r2, pc2 );

  circle_imp_print_2d ( r2, pc2, "  The recovered implicit circle:" );

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test016 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST016 tests CIRCLE_IMP_POINTS_2D and POLYGON_AREA_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  int n;
  double *p;
  double pc[DIM_NUM] = { 5.0, -2.0 };
  double pi = 3.141592653589793;
  double r = 2.0;
  double result;

  cout << "\n";
  cout << "TEST016\n";
  cout << "  CIRCLE_IMP_POINTS_2D gets points on a circle;\n";
  cout << "  POLYGON_AREA_2D finds the area of a polygon.\n";

  circle_imp_print_2d ( r, pc, "  The implicit circle:" );

  cout << "\n";
  cout << "  The area = " << pi * r * r << "\n";

  n = 8;

  p = circle_imp_points_2d ( r, pc, n );

  r8mat_transpose_print ( DIM_NUM, n, p, "  Sample results:" );

  delete [] p;

  cout << "\n";
  cout << "  For any N, the sampled points define a polygon\n";
  cout << "  whose area approximates the circle area.\n";
  cout << "\n";
  cout << "  N      Area\n";
  cout << "\n";

  for ( n = 3; n <= 24; n++ )
  {
    p = circle_imp_points_2d ( r, pc, n );
    result = polygon_area_2d ( n, p );
    cout << "  " << setw(6) << n
         << "  " << setw(12) << result << "\n";

    delete [] p;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0165 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0165 tests CIRCLE_IMP_POINTS_3D.
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
{
# define DIM_NUM 3
# define N 12

  double nc[DIM_NUM] = { 1.0, 1.0, 1.0 };
  double *p;
  double pc[DIM_NUM] = { 5.0, -2.0, 1.0 };
  double r = 2.0;

  cout << "\n";
  cout << "TEST0165\n";
  cout << "  CIRCLE_IMP_POINTS_3D gets points on a circle in 3D;\n";

  circle_imp_print_3d ( r, pc, nc, "  The implicit circle:" );

  p = circle_imp_points_3d ( r, pc, nc, N );

  r8mat_transpose_print ( DIM_NUM, N, p, "  Points on the circle:" );

  delete [] p;

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test017 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST017 tests CIRCLE_IMP_POINTS_ARC_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 13
# define DIM_NUM 2

  double r = 2.0;
  double theta1;
  double theta2;
  double p[DIM_NUM*N];
  double pc[DIM_NUM] = { 5.0, -2.0 };
  double pi = 3.141592653589793;

  theta1 = pi / 2.0;
  theta2 = 3.0 * pi / 2.0;

  cout << "\n";
  cout << "TEST017\n";
  cout << "  CIRCLE_IMP_POINTS_ARC_2D returns points on a\n";
  cout << "  circular arc.\n";
  cout << "\n";
  cout << "  The circle will have center " << pc[0] << ", "<< pc[1] << "\n";
  cout << "  and radius R = " << r << "\n";
  cout << "\n";
  cout << "  The arc extends from THETA1 = " << theta1 << "\n";
  cout << "  to THETA2 = " << theta2 << "\n";

  circle_imp_points_arc_2d ( r, pc, theta1, theta2, N, p );

  r8mat_transpose_print ( DIM_NUM, N, p, "  Sample results:" );

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test018 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST018 tests CIRCLE_IMP_POINT_DIST_2D and CIRCLES_IMP_INT_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 5

  double d1;
  double d2;
  int int_num;
  double pint[DIM_NUM*2];
  double pc1[DIM_NUM] = { 0.0, 0.0 };
  double *pc2;
  double pc2_test[DIM_NUM*TEST_NUM] = {
    5.0,       5.0,
    7.0710678, 7.0710678,
    4.0,       0.0,
    6.0,       0.0,
    0.0,       0.0 };
  double r1 = 5.0;
  double r2;
  double r2_test[TEST_NUM] = { 0.5, 5.0, 3.0, 3.0, 5.0 };
  int test;

  cout << "\n";
  cout << "TEST018\n";
  cout << "  CIRCLE_IMP_POINT_DIST_2D finds the\n";
  cout << "    distance from a point to a circle.\n";
  cout << "  CIRCLES_IMP_INT_2D determines the intersections of\n";
  cout << "    two circles in 2D.\n";

  circle_imp_print_2d ( r1, pc1, "  The first circle:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    r2 = r2_test[test];
    pc2 = pc2_test + test*DIM_NUM;

    circle_imp_print_2d ( r2, pc2, "  The second circle:" );

    circles_imp_int_2d ( r1, pc1, r2, pc2, &int_num, pint );

    if ( int_num == 0 )
    {
      cout << "\n";
      cout << "  The circles do not intersect.\n";
    }
    else if ( int_num == 1 )
    {
      cout << "\n";
      cout << "  The circles intersect at one point:\n";
      cout << "\n";
      cout << "        P         Dist 1  Dist 2\n";
      cout << "\n";
      d1 = circle_imp_point_dist_2d ( r1, pc1, pint+0*DIM_NUM );
      d2 = circle_imp_point_dist_2d ( r2, pc2, pint+0*DIM_NUM );
      cout << "  " << setw(8) << pint[0+0*DIM_NUM]
           << "  " << setw(8) << pint[1+0*DIM_NUM]
           << "  " << setw(8) << d1
           << "  " << setw(8) << d2 << "\n";
    }
    else if ( int_num == 2 )
    {
      cout << "\n";
      cout << "  The circles intersect at two points:\n";
      cout << "\n";
      cout << "        P         Dist 1  Dist 2\n";
      cout << "\n";
      d1 = circle_imp_point_dist_2d ( r1, pc1, pint+0*DIM_NUM );
      d2 = circle_imp_point_dist_2d ( r2, pc2, pint+0*DIM_NUM );
      cout << "  " << setw(8) << pint[0+0*DIM_NUM]
           << "  " << setw(8) << pint[1+0*DIM_NUM]
           << "  " << setw(8) << d1
           << "  " << setw(8) << d2 << "\n";
      d1 = circle_imp_point_dist_2d ( r1, pc1, pint+1*DIM_NUM );
      d2 = circle_imp_point_dist_2d ( r2, pc2, pint+1*DIM_NUM );
      cout << "  " << setw(8) << pint[0+1*DIM_NUM]
           << "  " << setw(8) << pint[1+1*DIM_NUM]
           << "  " << setw(8) << d1
           << "  " << setw(8) << d2 << "\n";
    }
    else if ( int_num == 3 )
    {
      cout << "\n";
      cout << "  The circles coincide (infinite intersection).\n";
    }
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0183 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0183 tests CIRCLE_LLR2IMP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double d1;
  double d2;
  double d3;
  double d4;
  double p_hi =  10.0;
  double p_lo = -10.0;
  double *pc;
  double *p1;
  double *p2;
  double *q1;
  double *q2;
  double r;
  double r_hi;
  double r_lo;
  int seed;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "TEST0183\n";
  cout << "  CIRCLE_LLR2IMP_3D is given:\n";
  cout << "  a line through P1 and P2,\n";
  cout << "  a line through Q1 and Q2,\n";
  cout << "  and a radius R,\n";
  cout << "  and determines the centers C of 4 circles\n";
  cout << "  of the given radius, tangent to both lines.\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    p1 = r8vec_uniform_new ( DIM_NUM, p_lo, p_hi, seed );
    p2 = r8vec_uniform_new ( DIM_NUM, p_lo, p_hi, seed );
    q1 = r8vec_uniform_new ( DIM_NUM, p_lo, p_hi, seed );
    q2 = r8vec_uniform_new ( DIM_NUM, p_lo, p_hi, seed );

    r_lo = 1.0;
    r_hi = 5.0;
    r = r8_uniform ( r_lo, r_hi, seed );

    cout << "\n";
    cout << "  Radius R = " << r << "\n";

    cout << "  Point #P1: ( "
         << setw(12) << p1[0] << ","
         << setw(12) << p1[1] << ")\n";
    cout << "  Point #P2: ( "
         << setw(12) << p2[0] << ","
         << setw(12) << p2[1] << ")\n";
    cout << "  Point #Q1: ( "
         << setw(12) << q1[0] << ","
         << setw(12) << q1[1] << ")\n";
    cout << "  Point #Q2: ( "
         << setw(12) << q2[0] << ","
         << setw(12) << q2[1] << ")\n";

    pc = circle_llr2imp_2d ( p1, p2, q1, q2, r );

    cout << "  Center #1: ( "
         << setw(12) << pc[0+0*DIM_NUM] << ","
         << setw(12) << pc[1+0*DIM_NUM] << ")\n";
    cout << "  Center #2: ( "
         << setw(12) << pc[0+1*DIM_NUM] << ","
         << setw(12) << pc[1+1*DIM_NUM] << ")\n";
    cout << "  Center #3: ( "
         << setw(12) << pc[0+2*DIM_NUM] << ","
         << setw(12) << pc[1+2*DIM_NUM] << ")\n";
    cout << "  Center #4: ( "
         << setw(12) << pc[0+3*DIM_NUM] << ","
         << setw(12) << pc[1+3*DIM_NUM] << ")\n";
//
//  Check that the lines are the right distance from the center.
//
    d1 = line_exp_point_dist_2d ( p1, p2, pc+0*2 );
    d2 = line_exp_point_dist_2d ( p1, p2, pc+1*2 );
    d3 = line_exp_point_dist_2d ( p1, p2, pc+2*2 );
    d4 = line_exp_point_dist_2d ( p1, p2, pc+3*2 );

    cout << "  " << setw(12) << d1
         << "  " << setw(12) << d2
         << "  " << setw(12) << d3
         << "  " << setw(12) << d4 << "\n";

    d1 = line_exp_point_dist_2d ( q1, q2, pc+0*2 );
    d2 = line_exp_point_dist_2d ( q1, q2, pc+1*2 );
    d3 = line_exp_point_dist_2d ( q1, q2, pc+2*2 );
    d4 = line_exp_point_dist_2d ( q1, q2, pc+3*2 );

    cout << "  " << setw(12) << d1
         << "  " << setw(12) << d2
         << "  " << setw(12) << d3
         << "  " << setw(12) << d4 << "\n";

    delete [] p1;
    delete [] p2;
    delete [] pc;
    delete [] q1;
    delete [] q2;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0185 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0185 tests CIRCLE_PPPR2IMP_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double d11;
  double d12;
  double d21;
  double d22;
  int i;
  double normal[DIM_NUM];
  double p_hi =  10.0;
  double p_lo = -10.0;
  double pc[DIM_NUM*2];
  double *p1;
  double *p2;
  double *p3;
  double r;
  double r_hi;
  double r_lo;
  int seed;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "TEST0185\n";
  cout << "  CIRCLE_PPPR2IMP_3D is given 3D points P1, P2, P3,\n";
  cout << "  and a radius R,\n";
  cout << "  and determines the centers C of two circles\n";
  cout << "  of the given radius, passing through P1 and P2\n";
  cout << "  and lying in the plane of P1, P2 and P3.\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    p1 = r8vec_uniform_new ( DIM_NUM, p_lo, p_hi, seed );
    p2 = r8vec_uniform_new ( DIM_NUM, p_lo, p_hi, seed );
    p3 = r8vec_uniform_new ( DIM_NUM, p_lo, p_hi, seed );

    r_lo = r8vec_distance ( DIM_NUM, p1, p2 );
    r_hi = r_lo + 5.0;
    r = r8_uniform ( r_lo, r_hi, seed );

    cout << "\n";
    cout << "  Radius R = " << r << "\n";

    cout << "  Point #1: ( "
         << setw(12) << p1[0] << ","
         << setw(12) << p1[1] << ")\n";
    cout << "  Point #2: ( "
         << setw(12) << p2[0] << ","
         << setw(12) << p2[1] << ")\n";
    cout << "  Point #3: ( "
         << setw(12) << p3[0] << ","
         << setw(12) << p3[1] << ")\n";

    circle_pppr2imp_3d ( p1, p2, p3, r, pc, normal );

    cout << "  Center #1: ( "
         << setw(12) << pc[0+0*DIM_NUM] << ","
         << setw(12) << pc[1+0*DIM_NUM] << ")\n";
    cout << "  Center #2: ( "
         << setw(12) << pc[0+1*DIM_NUM] << ","
         << setw(12) << pc[1+1*DIM_NUM] << ")\n";
//
//  Check that the points are the right distance from the center.
//
    d11 = r8vec_distance ( DIM_NUM, p1, pc+0*DIM_NUM );
    d21 = r8vec_distance ( DIM_NUM, p2, pc+0*DIM_NUM );
    d12 = r8vec_distance ( DIM_NUM, p1, pc+1*DIM_NUM );
    d22 = r8vec_distance ( DIM_NUM, p2, pc+1*DIM_NUM );

    cout << "  " << setw(12) << d11
         << "  " << setw(12) << d21
         << "  " << setw(12) << d12
         << "  " << setw(12) << d22 << "\n";
//
//  Check that the radial vector to the point is perpendicular to NORMAL.
//
    d11 = 0.0;
    d21 = 0.0;
    d12 = 0.0;
    d22 = 0.0;
    for ( i = 0; i < DIM_NUM; i++ )
    {
      d11 = d11 + normal[i] * ( p1[i] - pc[i+0*DIM_NUM] );
      d21 = d21 + normal[i] * ( p2[i] - pc[i+0*DIM_NUM] );
      d12 = d12 + normal[i] * ( p1[i] - pc[i+1*DIM_NUM] );
      d22 = d22 + normal[i] * ( p2[i] - pc[i+1*DIM_NUM] );
    }

    cout << "  " << setw(12) << d11
         << "  " << setw(12) << d21
         << "  " << setw(12) << d12
         << "  " << setw(12) << d22 << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test019 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST019 tests CIRCLE_PPR2IMP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 November 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double d11;
  double d12;
  double d21;
  double d22;
  double p_hi =  10.0;
  double p_lo = -10.0;
  double *pc;
  double *p1;
  double *p2;
  double r;
  double r_hi;
  double r_lo;
  int seed;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "TEST019\n";
  cout << "  CIRCLE_PPR2IMP_2D is given points P1 and P2,\n";
  cout << "  and a radius R,\n";
  cout << "  and determines the centers C of two circles\n";
  cout << "  of the given radius, passing through P1 and P2.\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    p1 = r8vec_uniform_new ( DIM_NUM, p_lo, p_hi, seed );
    p2 = r8vec_uniform_new ( DIM_NUM, p_lo, p_hi, seed );

    r_lo = r8vec_distance ( DIM_NUM, p1, p2 );
    r_hi = r_lo + 5.0;
    r = r8_uniform ( r_lo, r_hi, seed );

    cout << "\n";
    cout << "  Radius R = " << r << "\n";

    cout << "  Point #1: ( "
         << setw(12) << p1[0] << ","
         << setw(12) << p1[1] << ")\n";
    cout << "  Point #2: ( "
         << setw(12) << p2[0] << ","
         << setw(12) << p2[1] << ")\n";

    pc = circle_ppr2imp_2d ( p1, p2, r );

    cout << "  Center #1: ( "
         << setw(12) << pc[0+0*DIM_NUM] << ","
         << setw(12) << pc[1+0*DIM_NUM] << ")\n";
    cout << "  Center #2: ( "
         << setw(12) << pc[0+1*DIM_NUM] << ","
         << setw(12) << pc[1+1*DIM_NUM] << ")\n";

    d11 = r8vec_distance ( DIM_NUM, p1, pc+0*DIM_NUM );
    d21 = r8vec_distance ( DIM_NUM, p2, pc+0*DIM_NUM );
    d12 = r8vec_distance ( DIM_NUM, p1, pc+1*DIM_NUM );
    d22 = r8vec_distance ( DIM_NUM, p2, pc+1*DIM_NUM );

    cout << "  " << setw(12) << d11
         << "  " << setw(12) << d21
         << "  " << setw(12) << d12
         << "  " << setw(12) << d22 << "\n";

    delete [] p1;
    delete [] p2;
    delete [] pc;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test020 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST020 tests CUBE_SIZE_3D and CUBE_SHAPE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int edge_num;
  int face_num;
  int *face_order;
  int face_order_max;
  int *face_point;
  int point_num;
  double *point_coord;

  cout << "\n";
  cout << "TEST020\n";
  cout << "  For the cube,\n";
  cout << "  CUBE_SIZE_3D returns dimension information;\n";
  cout << "  CUBE_SHAPE_3D returns face and order information.\n";
  cout << "  SHAPE_PRINT_3D prints this information.\n";
//
//  Get the sizes.
//
  cube_size_3d ( &point_num, &edge_num, &face_num, &face_order_max );

  cout << "\n";
  cout << "    Number of vertices: " << point_num << "\n";
  cout << "    Number of edges   : " << edge_num << "\n";
  cout << "    Number of faces   : " << face_num << "\n";
  cout << "    Maximum face order: " << face_order_max << "\n";
//
//  Make room for the data.
//
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];
  point_coord = new double[DIM_NUM*point_num];
//
//  Get the data.
//
  cube_shape_3d ( point_num, face_num, face_order_max, point_coord,
    face_order, face_point );
//
//  Print the data.
//
  shape_print_3d ( point_num, face_num, face_order_max,
    point_coord, face_order, face_point );

  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0201 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0201 tests CYLINDER_POINT_DIST_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 6

  double dist;
  double dist_test[TEST_NUM] = {
    3.0, 0.5, 5.0, 8.0, 1.0, 0.25 };
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
      4.0,    0.5,  0.0,
     -0.5,   -1.0,  0.0,
      4.0,    6.0,  0.0,
      0.75, -10.0,  0.0,
      0.0,    0.0,  0.0,
      0.25,   1.75, 0.0 };
  double p1[DIM_NUM] = { 0.0, -2.0, 0.0 };
  double p2[DIM_NUM] = { 0.0,  2.0, 0.0 };
  double r = 1.0;
  int test;

  cout << "\n";
  cout << "TEST0201\n";
  cout << "  CYLINDER_POINT_DIST_3D computes the distance\n";
  cout << "  to a cylinder.\n";

  cout << "\n";
  cout << "  Radius R = " << r << "\n";
  cout << "  Center of bottom disk ="
       << "  " << setw(10) << p1[0]
       << "  " << setw(10) << p1[1]
       << "  " << setw(10) << p1[2] << "\n";
  cout << "  Center of top disk =   "
       << "  " << setw(10) << p2[0]
       << "  " << setw(10) << p2[1]
       << "  " << setw(10) << p2[2] << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    cout << "\n";
    cout << "  P ="
         << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << p[2] << "\n";

    dist = cylinder_point_dist_3d ( p1, p2, r, p );

    cout << "  Distance (computed) ="
         << "  " << setw(10) << dist << "\n";
    cout << "  Distance (exact)    ="
         << "  " << setw(10) << dist_test[test] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test02015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02015 tests CYLINDER_POINT_DIST_SIGNED_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 6

  double dist;
  double dist_test[TEST_NUM] = {
    3.0, -0.5, 5.0, 8.0, -1.0, -0.25 };
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
      4.0,    0.5,  0.0,
     -0.5,   -1.0,  0.0,
      4.0,    6.0,  0.0,
      0.75, -10.0,  0.0,
      0.0,    0.0,  0.0,
      0.25,   1.75, 0.0 };
  double p1[DIM_NUM] = { 0.0, -2.0, 0.0 };
  double p2[DIM_NUM] = { 0.0,  2.0, 0.0 };
  double r = 1.0;
  int test;

  cout << "\n";
  cout << "TEST02015\n";
  cout << "  CYLINDER_POINT_DIST_SIGNED_3D computes the signed\n";
  cout << "  distance to a cylinder.\n";

  cout << "\n";
  cout << "  Radius R = " << r << "\n";
  cout << "  Center of bottom disk ="
       << "  " << setw(10) << p1[0]
       << "  " << setw(10) << p1[1]
       << "  " << setw(10) << p1[2] << "\n";
  cout << "  Center of top disk =   "
       << "  " << setw(10) << p2[0]
       << "  " << setw(10) << p2[1]
       << "  " << setw(10) << p2[2] << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    cout << "\n";
    cout << "  P ="
         << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << p[2] << "\n";

    dist = cylinder_point_dist_signed_3d ( p1, p2, r, p );

    cout << "  Signed distance (computed) ="
         << "  " << setw(10) << dist << "\n";
    cout << "  Signed distance (exact)    ="
         << "  " << setw(10) << dist_test[test] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0202 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0202 tests CYLINDER_POINT_INSIDE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 6

  bool inside;
  bool inside_test[TEST_NUM] = {
    false, true, false, false, true, true };
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
      4.0,    0.5,  0.0,
     -0.5,   -1.0,  0.0,
      4.0,    6.0,  0.0,
      0.75, -10.0,  0.0,
      0.0,    0.0,  0.0,
      0.25,   1.75, 0.0 };
  double p1[DIM_NUM] = { 0.0, -2.0, 0.0 };
  double p2[DIM_NUM] = { 0.0,  2.0, 0.0 };
  double r = 1.0;
  int test;

  cout << "\n";
  cout << "TEST0202\n";
  cout << "  CYLINDER_POINT_INSIDE_3D determines if a point\n";
  cout << "  is inside a cylinder.\n";

  cout << "\n";
  cout << "  Radius R = " << r << "\n";
  cout << "  Center of bottom disk ="
       << "  " << setw(10) << p1[0]
       << "  " << setw(10) << p1[1]
       << "  " << setw(10) << p1[2] << "\n";
  cout << "  Center of top disk =   "
       << "  " << setw(10) << p2[0]
       << "  " << setw(10) << p2[1]
       << "  " << setw(10) << p2[2] << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    cout << "\n";
    cout << "  P ="
         << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << p[2] << "\n";

    inside = cylinder_point_inside_3d ( p1, p2, r, p );

    cout << "  INSIDE (computed) ="
         << "  " << setw(1) << inside << "\n";
    cout << "  INSIDE (exact)    ="
         << "  " << setw(1) << inside_test[test] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0203 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0203 tests CYLINDER_POINT_NEAR_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 6

  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
      4.0,    0.5,  0.0,
     -0.5,   -1.0,  0.0,
      4.0,    6.0,  0.0,
      0.75, -10.0,  0.0,
      0.0,    0.0,  0.0,
      0.25,   1.75, 0.0 };
  double p1[DIM_NUM] = { 0.0, -2.0, 0.0 };
  double p2[DIM_NUM] = { 0.0,  2.0, 0.0 };
  double *pn;
  double pn_test[DIM_NUM*TEST_NUM] = {
      1.0,   0.5,  0.0,
     -1.0,  -1.0,  0.0,
      1.0,   2.0,  0.0,
      0.75, -2.0,  0.0,
      1.0,   0.0,  0.0,
      0.25,  2.0,  0.0 };
  double r = 1.0;
  int test;

  cout << "\n";
  cout << "TEST0203\n";
  cout << "  CYLINDER_POINT_NEAR_3D computes the nearest point\n";
  cout << "  on a cylinder.\n";

  cout << "\n";
  cout << "  Radius R = " << r << "\n";
  cout << "  Center of bottom disk ="
       << "  " << setw(10) << p1[0]
       << "  " << setw(10) << p1[1]
       << "  " << setw(10) << p1[2] << "\n";
  cout << "  Center of top disk =   "
       << "  " << setw(10) << p2[0]
       << "  " << setw(10) << p2[1]
       << "  " << setw(10) << p2[2] << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    cout << "\n";
    cout << "  P ="
         << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << p[2] << "\n";

    pn = cylinder_point_near_3d ( p1, p2, r, p );

    cout << "  PN (computed) ="
         << "  " << setw(10) << pn[0]
         << "  " << setw(10) << pn[1]
         << "  " << setw(10) << pn[2] << "\n";
    cout << "  PN (exact)    ="
         << "  " << setw(10) << pn_test[0+test*DIM_NUM]
         << "  " << setw(10) << pn_test[1+test*DIM_NUM]
         << "  " << setw(10) << pn_test[2+test*DIM_NUM] << "\n";
    delete [] pn;
  }

  cout << "\n";
  cout << "  (Note that case 5 is ambiguous.  The set of nearest\n";
  cout << "  points forms a circle, any of which will do.)\n";

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test02035 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02035 tests CYLINDER_SAMPLE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define N 20

  double *p;
  double p1[DIM_NUM] = { 0.0, -2.0, 0.0 };
  double p2[DIM_NUM] = { 0.0,  2.0, 0.0 };
  double r = 1.0;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST02035\n";
  cout << "  CYLINDER_SAMPLE_3D samples points in a cylinder.\n";

  cout << "\n";
  cout << "  Radius R = " << r << "\n";
  cout << "  Center of bottom disk ="
       << "  " << setw(10) << p1[0]
       << "  " << setw(10) << p1[1]
       << "  " << setw(10) << p1[2] << "\n";
  cout << "  Center of top disk =   "
       << "  " << setw(10) << p2[0]
       << "  " << setw(10) << p2[1]
       << "  " << setw(10) << p2[2] << "\n";

  p = cylinder_sample_3d ( p1, p2, r, N, seed );

  r8mat_transpose_print ( DIM_NUM, N, p, "  Sample points:" );

  delete [] p;

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test0204 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0204 tests CYLINDER_VOLUME_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double p1[DIM_NUM] = { 1.0, 2.0, 3.0 };
  double p2[DIM_NUM] = { 5.0, 6.0, 5.0 };
  double pi = 3.141592653589793;
  double r = 5.0;
  double volume;

  cout << "\n";
  cout << "TEST0204\n";
  cout << "  CYLINDER_VOLUME_3D computes the volume of a cylinder.\n";

  cout << "\n";
  cout << "  Radius R = " << r << "\n";
  cout << "  Center of bottom disk ="
       << "  " << setw(10) << p1[0]
       << "  " << setw(10) << p1[1]
       << "  " << setw(10) << p1[2] << "\n";
  cout << "  Center of top disk =   "
       << "  " << setw(10) << p2[0]
       << "  " << setw(10) << p2[1]
       << "  " << setw(10) << p2[2] << "\n";

  volume = cylinder_volume_3d ( p1, p2, r );

  cout << "\n";
  cout << "  Volume (computed) = " << volume << "\n";
  cout << "  Volume (exact)    = " << pi * 150.0 << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0205 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0205 tests DEGREES_TO_RADIANS and RADIANS_TO_DEGREES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double angle_deg;
  double angle_deg2;
  double angle_rad;
  int i;

  cout << "\n";
  cout << "TEST0205\n";
  cout << "  DEGREES_TO_RADIANS converts an angle from degrees\n";
  cout << "    to radians;\n";
  cout << "  RADIANS_TO_DEGREES converts an angle from radians\n";
  cout << "    to degrees;\n";
  cout << "\n";
  cout << "  Degrees     Radians     Degrees\n";
  cout << "\n";

  for ( i = -2; i <= 14; i++ )
  {
    angle_deg = ( double ) ( 30 * i );
    angle_rad = degrees_to_radians ( angle_deg );
    angle_deg2 = radians_to_degrees ( angle_rad );
    cout << "  " << setw(10) << angle_deg
         << "  " << setw(10) << angle_rad
         << "  " << setw(10) << angle_deg2 << "\n";
  }

  return;
}
//****************************************************************************80

void test021 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST021 tests DIRECTION_PERT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int i;
  int itest;
  int seed;
  double sigma;
  double vbase[DIM_NUM] = { 1.0, 0.0, 0.0 };
  double *vran;

  cout << "\n";
  cout << "TEST021\n";
  cout << "  DIRECTION_PERT_3D perturbs a direction vector.\n";
  cout << "\n";

  seed = get_seed ( );

  cout << "\n";
  cout << "  We use SEED = " << seed << "\n";

  cout << "\n";
  cout << "  Base vector:\n";
  cout << "  " << vbase[0]
       << "  " << vbase[1]
       << "  " << vbase[2] << "\n";

  for ( itest = 0; itest < 3; itest++ )
  {
    if ( itest == 0 )
    {
      sigma = 0.99;
    }
    else if ( itest == 1 )
    {
      sigma = 0.5;
    }
    else
    {
      sigma = 0.1;
    }

    cout << "\n";
    cout << "  Using sigma = " << sigma << "\n";
    cout << "\n";

    for ( i = 0; i < 20; i++ )
    {
      vran = direction_pert_3d ( sigma, vbase, seed );
      cout << "  " << setw(10) << vran[0]
           << "  " << setw(10) << vran[1]
           << "  " << setw(10) << vran[2] << "\n";
      delete [] vran;
    }

  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test022 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST022 tests DIRECTION_UNIFORM_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  double *vran;

  cout << "\n";
  cout << "TEST022\n";
  cout << "  DIRECTION_UNIFORM_3D picks a random direction vector.\n";

  seed = get_seed ( );

  cout << "\n";
  cout << "  We use SEED = " << seed << "\n";

  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    vran = direction_uniform_3d ( seed );
    cout << "  " << setw(10) << vran[0]
         << "  " << setw(10) << vran[1]
         << "  " << setw(10) << vran[2] << "\n";
    delete [] vran;
  }

  return;
}
//****************************************************************************80

void test023 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST023 tests DIRECTION_UNIFORM_ND and GET_SEED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 4

  int i;
  int j;
  int seed = 123456789;
  double *vran;

  cout << "\n";
  cout << "TEST023\n";
  cout << "  DIRECTION_UNIFORM_ND picks a random direction vector.\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    vran = direction_uniform_nd ( DIM_NUM, seed );
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << "  " << setw(8) << vran[j];
    }
    cout << "\n";
    delete [] vran;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0232 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0232 tests DISK_POINT_DIST_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 5

  double axis[DIM_NUM] = { 0.0, 1.0, 1.0 };
  double dist;
  double dist_test[TEST_NUM] = { 2.0, 0.0, 0.0, 8.0, 10.0 };
  double *p;
  double pc[DIM_NUM] = { 0.0, 1.4142135, 1.4142135 };
  double p_test[DIM_NUM*TEST_NUM] = {
     0.0,  0.0,        0.0,
     0.0,  0.70710677, 2.1213202,
     2.0,  1.4142135,  1.4142135,
    10.0,  1.4142135,  1.4142135,
    10.0,  5.6568542,  5.6568542 };
  double r = 2.0;
  int test;

  cout << "\n";
  cout << "TEST0232\n";
  cout << "  DISK_POINT_DIST_3D finds the distance from\n";
  cout << "    a disk to a point in 3D.\n";

  cout << "\n";
  cout << "  Disk radius = " << r << "\n";
  r8vec_print ( DIM_NUM, pc, "  Disk center: " );
  r8vec_print ( DIM_NUM, axis, "  Disk axis: " );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    r8vec_print ( DIM_NUM, p, "  Point: " );

    dist = disk_point_dist_3d ( pc, r, axis, p );

    cout << "\n";
    cout << "  Distance = " << dist
         << "  Expected = " << dist_test[test] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0234 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0234 tests R8MAT_SOLVE_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double det;
  int i;
  int n = 2;
  int seed;
  int test;
  int test_num = 5;
  double *x;
  double *x2;

  cout << "\n";
  cout << "TEST0234\n";
  cout << "  R8MAT_SOLVE_2D solves 2D linear systems.\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r8mat_uniform_01_new ( n, n, seed );
    x = r8vec_uniform_01_new ( n, seed );
    b = r8mat_mv ( n, n, a, x );

    x2 = r8mat_solve_2d ( a, b, &det );

    cout << "\n";
    cout << "  Solution / Computed:\n";
    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(14) << x[i]
           << "  " << setw(14) << x2[i] << "\n";
    }

    delete [] a;
    delete [] b;
    delete [] x;
    delete [] x2;
  }

  return;
}
//****************************************************************************80

void test0235 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0235 tests DMS_TO_RADIANS and RADIANS_TO_DMS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int angle_deg;
  int angle_min;
  double angle_rad;
  double angle_rad2;
  int angle_sec;
  int i;
  double pi = 3.141592653589793;

  cout << "\n";
  cout << "TEST0235\n";
  cout << "  DMS_TO_RADIANS converts an angle from\n";
  cout << "    degrees/minutes/seconds to radians;\n";
  cout << "  RADIANS_TO_DEGREES converts an angle from radians\n";
  cout << "    to degrees/minutes/seconds;\n";
  cout << "\n";
  cout << "  Radians     DMS     Radians\n";
  wcout << "\n";

  for ( i = -2; i <= 15; i++ )
  {
    angle_rad = pi * ( double ) ( i ) / 7.0;

    radians_to_dms ( angle_rad, &angle_deg, &angle_min, &angle_sec );

    angle_rad2 = dms_to_radians ( angle_deg, angle_min, angle_sec );

    cout << "  " << setw(10) << angle_rad
         << "  " << setw(4)  << angle_deg
         << "  " << setw(3)  << angle_min
         << "  " << setw(3)  << angle_sec
         << "  " << setw(10) << angle_rad2 << "\n";
  }

  return;
}
//****************************************************************************80

void test0236 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0236 tests DODEC_SIZE_3D, DODEC_SHAPE_3D, SHAPE_PRINT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int edge_num;
  int face_num;
  int *face_order;
  int face_order_max;
  int *face_point;
  int point_num;
  double *point_coord;

  cout << "\n";
  cout << "TEST0236\n";
  cout << "  For the dodecahedron,\n";
  cout << "  DODEC_SIZE_3D returns dimension information;\n";
  cout << "  DODEC_SHAPE_3D returns face and order information.\n";
  cout << "  SHAPE_PRINT_3D prints this information.\n";
//
//  Get the sizes.
//
  dodec_size_3d ( &point_num, &edge_num, &face_num, &face_order_max );

  cout << "\n";
  cout << "    Number of vertices: " << point_num << "\n";
  cout << "    Number of edges   : " << edge_num << "\n";
  cout << "    Number of faces   : " << face_num << "\n";
  cout << "    Maximum face order: " << face_order_max << "\n";
//
//  Make room for the data.
//
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];
  point_coord = new double[DIM_NUM*point_num];
//
//  Get the data.
//
  dodec_shape_3d ( point_num, face_num, face_order_max, point_coord,
    face_order, face_point );
//
//  Print the data.
//
  shape_print_3d ( point_num, face_num, face_order_max,
    point_coord, face_order, face_point );

  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0238 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0238 tests DUAL_SIZE_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int edge_num1;
  int edge_num2;
  int *edge_point1;
  int face_num1;
  int face_num2;
  int *face_order1;
  int face_order_max1;
  int face_order_max2;
  int *face_point1;
  int point_num1;
  int point_num2;
  double *point_coord1;

  cout << "\n";
  cout << "TEST0238\n";
  cout << "  DUAL_SIZE_3D finds the sizes of the dual of a\n";
  cout << "  polyhedron;\n";
//
//  Get the CUBE shape.
//
  cube_size_3d ( &point_num1, &edge_num1, &face_num1, &face_order_max1 );

  cout << "\n";
  cout << "  The cube:\n";
  cout << "    Number of vertices: " << point_num1 << "\n";
  cout << "    Number of edges   : " << edge_num1 << "\n";
  cout << "    Number of faces   : " << face_num1 << "\n";
  cout << "    Maximum face order: " << face_order_max1 << "\n";

  face_order1 = new int[face_num1];
  face_point1 = new int[face_order_max1*face_num1];
  point_coord1 = new double[DIM_NUM*point_num1];

  cube_shape_3d ( point_num1, face_num1, face_order_max1, point_coord1,
    face_order1, face_point1 );

  dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, 
    point_coord1, face_order1, face_point1, &point_num2, &edge_num2,
    &face_num2, &face_order_max2 );

  cout << "\n";
  cout << "  The dual of the cube:\n";
  cout << "    Number of vertices: " << point_num2 << "\n";
  cout << "    Number of edges   : " << edge_num2 << "\n";
  cout << "    Number of faces   : " << face_num2 << "\n";
  cout << "    Maximum face order: " << face_order_max2 << "\n";

  delete [] face_order1;
  delete [] face_point1;
  delete [] point_coord1;
//
//  Get the DODECAHEDRON shape.
//
  dodec_size_3d ( &point_num1, &edge_num1, &face_num1, &face_order_max1 );

  cout << "\n";
  cout << "  The dodecahedron:\n";
  cout << "    Number of vertices: " << point_num1 << "\n";
  cout << "    Number of edges   : " << edge_num1 << "\n";
  cout << "    Number of faces   : " << face_num1 << "\n";
  cout << "    Maximum face order: " << face_order_max1 << "\n";

  face_order1 = new int[face_num1];
  face_point1 = new int[face_order_max1*face_num1];
  point_coord1 = new double[DIM_NUM*point_num1];

  dodec_shape_3d ( point_num1, face_num1, face_order_max1, point_coord1,
    face_order1, face_point1 );

  dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, 
    point_coord1, face_order1, face_point1, &point_num2, &edge_num2,
    &face_num2, &face_order_max2 );

  cout << "\n";
  cout << "  The dual of the dodecahedron:\n";
  cout << "    Number of vertices: " << point_num2 << "\n";
  cout << "    Number of edges   : " << edge_num2 << "\n";
  cout << "    Number of faces   : " << face_num2 << "\n";
  cout << "    Maximum face order: " << face_order_max2 << "\n";

  delete [] face_order1;
  delete [] face_point1;
  delete [] point_coord1;
//
//  Get the ICOSAHEDRON shape.
//
  icos_size ( &point_num1, &edge_num1, &face_num1, &face_order_max1 );

  cout << "\n";
  cout << "  The icosahedron:\n";
  cout << "    Number of vertices: " << point_num1 << "\n";
  cout << "    Number of edges   : " << edge_num1 << "\n";
  cout << "    Number of faces   : " << face_num1 << "\n";
  cout << "    Maximum face order: " << face_order_max1 << "\n";

  edge_point1 = new int[2*edge_num1];
  face_order1 = new int[face_num1];
  face_point1 = new int[face_order_max1*face_num1];
  point_coord1 = new double[DIM_NUM*point_num1];

  icos_shape ( point_num1, edge_num1, face_num1, face_order_max1, 
    point_coord1, edge_point1, face_order1, face_point1 );

  dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, 
    point_coord1, face_order1, face_point1, &point_num2, &edge_num2,
    &face_num2, &face_order_max2 );

  cout << "\n";
  cout << "  The dual of the icosahedron:\n";
  cout << "    Number of vertices: " << point_num2 << "\n";
  cout << "    Number of edges   : " << edge_num2 << "\n";
  cout << "    Number of faces   : " << face_num2 << "\n";
  cout << "    Maximum face order: " << face_order_max2 << "\n";

  delete [] edge_point1;
  delete [] face_order1;
  delete [] face_point1;
  delete [] point_coord1;
//
//  Get the octahedron shape.
//
  octahedron_size_3d ( &point_num1, &edge_num1, &face_num1, &face_order_max1 );

  cout << "\n";
  cout << "  The octahedron:\n";
  cout << "    Number of vertices: " << point_num1 << "\n";
  cout << "    Number of edges   : " << edge_num1 << "\n";
  cout << "    Number of faces   : " << face_num1 << "\n";
  cout << "    Maximum face order: " << face_order_max1 << "\n";

  face_order1 = new int[face_num1];
  face_point1 = new int[face_order_max1*face_num1];
  point_coord1 = new double[DIM_NUM*point_num1];

  octahedron_shape_3d ( point_num1, face_num1, face_order_max1,
    point_coord1, face_order1, face_point1 );

  dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, 
    point_coord1, face_order1, face_point1, &point_num2, &edge_num2,
    &face_num2, &face_order_max2 );

  cout << "\n";
  cout << "  The dual of the octahedron:\n";
  cout << "    Number of vertices: " << point_num2 << "\n";
  cout << "    Number of edges   : " << edge_num2 << "\n";
  cout << "    Number of faces   : " << face_num2 << "\n";
  cout << "    Maximum face order: " << face_order_max2 << "\n";

  delete [] face_order1;
  delete [] face_point1;
  delete [] point_coord1;
//
//  Get the soccer ball shape.
//
  soccer_size_3d ( &point_num1, &edge_num1, &face_num1, &face_order_max1 );

  cout << "\n";
  cout << "  The soccer ball:\n";
  cout << "    Number of vertices: " << point_num1 << "\n";
  cout << "    Number of edges   : " << edge_num1 << "\n";
  cout << "    Number of faces   : " << face_num1 << "\n";
  cout << "    Maximum face order: " << face_order_max1 << "\n";

  face_order1 = new int[face_num1];
  face_point1 = new int[face_order_max1*face_num1];
  point_coord1 = new double[DIM_NUM*point_num1];

  soccer_shape_3d ( point_num1, face_num1, face_order_max1,
    point_coord1, face_order1, face_point1 );

  dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, 
    point_coord1, face_order1, face_point1, &point_num2, &edge_num2,
    &face_num2, &face_order_max2 );

  cout << "\n";
  cout << "  The dual of the soccer ball:\n";
  cout << "    Number of vertices: " << point_num2 << "\n";
  cout << "    Number of edges   : " << edge_num2 << "\n";
  cout << "    Number of faces   : " << face_num2 << "\n";
  cout << "    Maximum face order: " << face_order_max2 << "\n";

  delete [] face_order1;
  delete [] face_point1;
  delete [] point_coord1;
//
//  Get the tetrahedron shape.
//
  tetrahedron_size_3d ( &point_num1, &edge_num1, &face_num1, &face_order_max1 );

  cout << "\n";
  cout << "  The tetrahedron:\n";
  cout << "    Number of vertices: " << point_num1 << "\n";
  cout << "    Number of edges   : " << edge_num1 << "\n";
  cout << "    Number of faces   : " << face_num1 << "\n";
  cout << "    Maximum face order: " << face_order_max1 << "\n";

  face_order1 = new int[face_num1];
  face_point1 = new int[face_order_max1*face_num1];
  point_coord1 = new double[DIM_NUM*point_num1];

  tetrahedron_shape_3d ( point_num1, face_num1, face_order_max1,
    point_coord1, face_order1, face_point1 );

  dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, 
    point_coord1, face_order1, face_point1, &point_num2, &edge_num2,
    &face_num2, &face_order_max2 );

  cout << "\n";
  cout << "  The dual of the tetrahedron:\n";
  cout << "    Number of vertices: " << point_num2 << "\n";
  cout << "    Number of edges   : " << edge_num2 << "\n";
  cout << "    Number of faces   : " << face_num2 << "\n";
  cout << "    Maximum face order: " << face_order_max2 << "\n";

  delete [] face_order1;
  delete [] face_point1;
  delete [] point_coord1;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test024 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST024 tests DUAL_SHAPE_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int edge_num1;
  int edge_num2;
  int face_num1;
  int face_num2;
  int *face_order1;
  int *face_order2;
  int face_order_max1;
  int face_order_max2;
  int *face_point1;
  int *face_point2;
  int point_num1;
  int point_num2;
  double *point_coord1;
  double *point_coord2;

  cout << "\n";
  cout << "TEST024\n";
  cout << "  DUAL_SHAPE_3D finds the dual of a polyhedron.\n";
//
//  Get the dodecahedron shape.
//
  cout << "\n";
  cout << "  The dodecahedron:\n";

  dodec_size_3d ( &point_num1, &edge_num1, &face_num1, &face_order_max1 );

  cout << "\n";
  cout << "    Number of vertices: " << point_num1 << "\n";
  cout << "    Number of edges   : " << edge_num1 << "\n";
  cout << "    Number of faces   : " << face_num1 << "\n";
  cout << "    Maximum face order: " << face_order_max1 << "\n";

  face_order1 = new int[face_num1];
  face_point1 = new int [face_order_max1*face_num1];
  point_coord1 = new double[DIM_NUM*point_num1];

  dodec_shape_3d ( point_num1, face_num1, face_order_max1, point_coord1,
    face_order1, face_point1 );

  shape_print_3d ( point_num1, face_num1, face_order_max1,
    point_coord1, face_order1, face_point1 );
//
//  Get the dual.
//
  cout << "\n";
  cout << "  The dual of the dodecahedron:\n";

  dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, 
    point_coord1, face_order1, face_point1, &point_num2, &edge_num2,
    &face_num2, &face_order_max2 );

  cout << "\n";
  cout << "    Number of vertices: " << point_num2 << "\n";
  cout << "    Number of edges   : " << edge_num2 << "\n";
  cout << "    Number of faces   : " << face_num2 << "\n";
  cout << "    Maximum face order: " << face_order_max2 << "\n";

  face_order2 = new int[face_num2];
  face_point2 = new int[face_order_max2*face_num2];
  point_coord2 = new double[DIM_NUM*point_num2];

  dual_shape_3d ( point_num1, face_num1, face_order_max1, point_coord1,
    face_order1, face_point1, point_num2, face_num2, face_order_max2,
    point_coord2, face_order2, face_point2 );

  shape_print_3d ( point_num2, face_num2, face_order_max2,
    point_coord2, face_order2, face_point2 );

  delete [] face_order1;
  delete [] face_order2;
  delete [] face_point1;
  delete [] face_point2;
  delete [] point_coord1;
  delete [] point_coord2;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0243 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0243 tests R8VEC_ANY_NORMAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 10
# define TEST_NUM 5

  int seed;
  int test;
  double *v1;
  double v1_norm;
  double v1v2_dot;
  double *v2;
  double v2_norm;

  cout << "\n";
  cout << "TEST0243\n";
  cout << "  R8VEC_ANY_NORMAL computes a vector V2 that is normal\n";
  cout << "  to a given vector V1.\n";
  cout << "\n";
  cout << "    Test    ||V1||      ||V2||        V1.V2\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 0; test < TEST_NUM; test++ )
  {
    v1 = r8vec_uniform_01_new ( DIM_NUM, seed );
    v1_norm = r8vec_norm ( DIM_NUM, v1 );
    v2 = r8vec_any_normal ( DIM_NUM, v1 );
    v2_norm = r8vec_norm ( DIM_NUM, v2 );
    v1v2_dot = r8vec_dot_product ( DIM_NUM, v1, v2 );
    cout << "  " << setw(6)  << test
         << "  " << setw(10) << v1_norm
         << "  " << setw(10) << v2_norm
         << "  " << setw(10) << v1v2_dot << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0245 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0245 tests R8VEC_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 1000

  int i;
  int j;
  int n;
  int seed = 123456789;
  double *x;
  double x_max;
  double x_mean;
  double x_min;
  double x_var;

  cout << "\n";
  cout << "TEST0245\n";
  cout << "  R8VEC_NORMAL_01 computes a vector of normally\n";
  cout << "  distributed random numbers.\n";
  cout << "  Using initial random number seed = " << seed << "\n";
//
//  Test 1:
//  Simply call 5 times for 1 value, and print.
//
  cout << "\n";
  cout << "  Test #1: Call 5 times, 1 value each time.\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    n = 1;
    x = r8vec_normal_01_new ( n, seed );
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(6)  << i
           << "  " << setw(10) << x[j] << "\n";
    }
    delete [] x;
  }
//
//  Test 2:
//  Restore the random number seed, and repeat.
//
  cout << "\n";
  cout << "  Test #2: Restore the random number seed.\n";
  cout << "  call 5 times, 1 value each time.\n";
  cout << "  The results should be identical.\n";
  cout << "\n";

  n = -1;
  x = r8vec_normal_01_new ( n, seed );
  delete [] x;

  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    n = 1;
    x = r8vec_normal_01_new ( n, seed );
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(6)  << i
           << "  " << setw(10) << x[j] << "\n";
    }
    delete [] x;
  }
//
//  Test 3:
//  Restore the random number seed, compute all 5 values at once.
//
  cout << "\n";
  cout << "  Test #3: Restore the random number seed.\n";
  cout << "  call 1 time for 5 values.\n";
  cout << "  The results should be identical.\n";
  cout << "\n";

  n = -1;
  x = r8vec_normal_01_new ( n, seed );
  delete [] x;

  seed = 123456789;

  n = 5;
  x = r8vec_normal_01_new ( n, seed );

  i = 0;
  for ( j = 0; j < n; j++ )
  {
    i = i + 1;
    cout << "  " << setw(6)  << i
         << "  " << setw(10) << x[j] << "\n";
  }
  delete [] x;
//
//  Test 4:
//  Restore the random number seed, compute all 5 values at once.
//
  cout << "\n";
  cout << "  Test #4: Restore the random number seed.\n";
  cout << "  call for 2, 1, and 2 values.\n";
  cout << "  The results should be identical.\n";
  cout << "\n";

  n = -1;
  x = r8vec_normal_01_new ( n, seed );
  delete [] x;

  seed = 123456789;

  n = 2;
  x = r8vec_normal_01_new ( n, seed );
  i = 0;

  for ( j = 0; j < n; j++ )
  {
    i = i + 1;
    cout << "  " << setw(6)  << i
         << "  " << setw(10) << x[j] << "\n";
  }
  delete [] x;

  n = 1;
  x = r8vec_normal_01_new ( n, seed );

  for ( j = 0; j < n; j++ )
  {
    i = i + 1;
    cout << "  " << setw(6)  << i
         << "  " << setw(10) << x[j] << "\n";
  }
  delete [] x;

  n = 2;
  x = r8vec_normal_01_new ( n, seed );

  for ( j = 0; j < n; j++ )
  {
    i = i + 1;
    cout << "  " << setw(6)  << i
         << "  " << setw(10) << x[j] << "\n";
  }
  delete [] x;
//
//  Test 5:
//  Determine the minimum, maximum, mean and variance.
//
  n = N_MAX;
  x = r8vec_normal_01_new ( n, seed );

  x_min = r8vec_min ( n, x );
  x_max = r8vec_max ( n, x );
  x_mean = r8vec_mean ( n, x );
  x_var = r8vec_variance ( n, x );

  cout << "\n";
  cout << "  Test #5:\n";
  cout << "  Number of samples was " << n       << "\n";
  cout << "  Minimum value was     " << x_min   << "\n";
  cout << "  Maximum value was     " << x_max   << "\n";
  cout << "  Average value was     " << x_mean  << "\n";
  cout << "  Variance was          " << x_var   << "\n";
  cout << "  Expected average      " << 0.0 << "\n";
  cout << "  Expected variance     " << 1.0 << "\n";

  delete [] x;

  return;
# undef N_MAX
}
//****************************************************************************80

void test025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025 tests ELLIPSE_POINT_DIST_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double dist;
  int i;
  int n = 10;
  double p[DIM_NUM];
  double r1 = 3.0;
  double r2 = 2.0;

  cout << "\n";
  cout << "TEST025:\n";
  cout << "  ELLIPSE_POINT_DIST_2D is given a point P, and\n";
  cout << "  finds the distance to an ellipse in 2D.\n";
  cout << "\n";
  cout << "  The ellipse is (X/R1)^2 + (Y/R2)^2 = 1\n";
  cout << "\n";
  cout << "  R1 = " << r1 << "\n";
  cout << "  R2 = " << r2 << "\n";
  cout << "\n";
  cout << "           P            DIST\n";
  cout << "\n";

  for ( i = -3; i <= n + 3; i++ )
  {
    p[0] = ( double ( n - i ) * 0.0
           + double (     i ) * 4.0 )
           / double ( n     );

    p[1] = ( double ( n - i ) * 3.0
           + double (     i ) * 0.0 )
           / double ( n     );

    dist = ellipse_point_dist_2d ( r1, r2, p );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  " << setw(8) << dist << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0255 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0255 tests ELLIPSE_POINT_NEAR_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  int i;
  int n = 10;
  double p[DIM_NUM];
  double *pn;
  double r1 = 3.0;
  double r2 = 2.0;

  cout << "\n";
  cout << "TEST0255:\n";
  cout << "  ELLIPSE_POINT_NEAR_2D is given a point P, and\n";
  cout << "  finds the nearest point on an ellipse in 2D.\n";
  cout << "\n";
  cout << "  The ellipse is (X/R1)^2 + (Y/R2)^2 = 1\n";
  cout << "\n";
  cout << "  R1 = " << r1 << "\n";
  cout << "  R2 = " << r2 << "\n";
  cout << "\n";
  cout << "           P            PN\n";
  cout << "\n";

  for ( i = -3; i <= n + 3; i++ )
  {
    p[0] = ( double ( n - i ) * 0.0
           + double (     i ) * 4.0 )
           / double ( n     );

    p[1] = ( double ( n - i ) * 3.0
           + double (     i ) * 0.0 )
           / double ( n     );

    pn = ellipse_point_near_2d ( r1, r2, p );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  " << setw(8) << pn[0]
         << "  " << setw(8) << pn[1] << "\n";

    delete [] pn;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test026 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026 tests ELLIPSE_POINTS_2D, ELLIPSE_AREA_2D, POLYGON_AREA_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N_MAX 24

  double area;
  int n;
  double pc[DIM_NUM] = { 5.0, -2.0 };
  double pi = 3.141592653589793;
  double psi = 3.141592653589793 / 6.0;
  double r1 = 3.0;
  double r2 = 1.0;
  double result;
  double *v;

  cout << "\n";
  cout << "TEST026\n";
  cout << "  ELLIPSE_POINTS_2D returns points on an ellipse;\n";
  cout << "  ELLIPSE_AREA_2D returns the area of an ellipse;\n";
  cout << "  POLYGON_AREA_2D finds the area of a polygon.\n";

  r8vec_print ( DIM_NUM, pc, "  Ellipse center:" );

  cout << "\n";
  cout << "  radii R1 = " << r1 << " R2 = " << r2 << "\n";
  cout << "  and angle PSI = " << psi << "\n";

  area = ellipse_area_2d ( r1, r2 );

  cout << "  and area = " << area << "\n";

  n = 16;
  v = new double[DIM_NUM*n];

  ellipse_points_2d ( pc, r1, r2, psi, n, v );

  r8mat_transpose_print ( DIM_NUM, n, v, "  Sample points:" );

  delete [] v;

  cout << "\n";
  cout << "  For any N, the sampled points define a polygon\n";
  cout << "  whose area approximates the ellipse area.\n";
  cout << "\n";
  cout << "       N         Area\n";
  cout << "\n";

  for ( n = 3; n <= N_MAX; n++ )
  {
    v = new double[DIM_NUM*n];
    ellipse_points_2d ( pc, r1, r2, psi, n, v );
    result = polygon_area_2d ( n, v );
    cout << "  " << setw(6) << n
         << "  " << setw(12) << result << "\n";
    delete [] v;
  }

  return;
# undef DIM_NUM
# undef N_MAX
}
//****************************************************************************80

void test027 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST027 tests ELLIPSE_POINTS_ARC_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 13

  double p[DIM_NUM*N];
  double pc[DIM_NUM] = { 5.0, -2.0 };
  double pi = 3.141592653589793;
  double psi = 3.141592653589793 / 6.0;
  double r1 = 3.0;
  double r2 = 1.0;
  double theta1 = 3.141592653589793 / 2.0;
  double theta2 = 2.0 * 3.141592653589793;

  cout << "\n";
  cout << "TEST027\n";
  cout << "  ELLIPSE_POINTS_ARC_2D returns points on an\n";
  cout << "    elliptical arc.\n";
  cout << "\n";
  cout << "  The ellipse has center " << pc[0] << "  " << pc[1] << "\n";
  cout << "  radii R1 = " << r1 << " R2 = " << r2 << "\n";
  cout << "  and angle PSI = " << psi << "\n";
  cout << "\n";
  cout << "  The arc extends from THETA1 = " << theta1 << "\n";
  cout << "  to THETA2 = " << theta2 << "\n";

  ellipse_points_arc_2d ( pc, r1, r2, psi, theta1, theta2, N, p );

  r8mat_transpose_print ( DIM_NUM, N, p, "  Sample points:" );

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test028 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST028 tests HALFPLANE_CONTAINS_POINT_2D
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  bool expected[TEST_NUM] = { true, false, true, false };
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
    1.0,   1.0,
    1.0,  -1.0,
   -1.0,   1.0,
    2.0, 200.0 };
  double *p1;
  double p1_test[DIM_NUM*TEST_NUM] = {
    0.0,   0.0,
    0.0,   0.0,
   -5.0,  -5.0,
    3.0, 150.0 };
  double *p2;
  double p2_test[DIM_NUM*TEST_NUM] = {
    2.0,   0.0,
    2.0,   0.0,
   10.0,  10.0,
    1.0,  50.0 };
  bool temp;
  int test;

  cout << "\n";
  cout << "TEST028\n";
  cout << "  HALFPLANE_CONTAINS_POINT_2D determines whether a\n";
  cout << "  halfplane bounded by PA:PB contains the\n";
  cout << "  point P.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p1_test + test*DIM_NUM;
    p2 = p2_test + test*DIM_NUM;
    p = p_test + test*DIM_NUM;

    temp = halfplane_contains_point_2d ( p1, p2, p );

    cout << "\n";
    cout << "  P1 = " << setw(12) << p1[0]
         << "  "      << setw(12) << p1[1] << "\n";
    cout << "  P2 = " << setw(12) << p2[0]
         << "  "      << setw(12) << p2[1] << "\n";
    cout << "  P = "  << setw(12) << p[0]
         << "  "      << setw(12) << p[1] << "\n";
    cout << "  Contains? = " << temp
         << "  Correct = " << expected[test] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test029 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST029 tests HALFSPACE_IMP_TRIANGLE_INT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double a;
  double b;
  double c;
  double d;
  int int_num;
  int j;
  int test;
  int test_num = 6;
  double p[DIM_NUM*4];
  double t[DIM_NUM*3];

  cout << "\n";
  cout << "TEST029\n";
  cout << "  HALFSPACE_IMP_TRIANGLE_INT_3D finds\n";
  cout << "    intersection points of an implicit\n";
  cout << "    halfspace and a triangle.\n";

  a =   1.0;
  b = - 2.0;
  c = - 3.0;
  d =   6.0;

  cout << "\n";
  cout << "  The implicitly defined bounding plane\n";
  cout << "  has the form: A*X + B*Y + C*Z + D = 0.\n";
  cout << "  A,B,C,D = " << a << "  " << b << "  " << c << "  " << d << "\n";

  for ( test = 0; test < test_num; test++ )
    {

    if ( test == 0 )
    {
      t[0+0*3] =  0.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] = -1.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] = -2.0;
    }
    else if ( test == 1 )
    {
      t[0+0*3] = -6.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] = -1.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] = -2.0;
    }
    else if ( test == 2 )
    {
      t[0+0*3] =  0.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] =  3.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] =  2.0;
    }
    else if ( test == 3 )
    {
      t[0+0*3] = -6.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] =  4.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] =  3.0;
    }
    else if ( test == 4 )
    {
      t[0+0*3] = -8.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] = -1.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] = -2.0;
    }
    else if ( test == 5 )
    {
      t[0+0*3] =  0.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] =  4.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] =  4.0;
    }

    cout << "\n";
    cout << "  Case " << test << "\n";

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    int_num = halfspace_imp_triangle_int_3d ( a, b, c, d, t, p );

    cout << "\n";
    cout << "  Number of intersection points is " << int_num << "\n";
    cout << "\n";

    if ( 0 < int_num )
    {
      for ( j = 0; j < int_num; j++ )
      {
        cout << "  " << setw(4)  << j
             << "  " << setw(10) << p[0+j*DIM_NUM]
             << "  " << setw(10) << p[1+j*DIM_NUM]
             << "  " << setw(10) << p[2+j*DIM_NUM] << "\n";
      }
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test030 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST030 tests HALFSPACE_NORM_TRIANGLE_INT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int int_num;
  int j;
  double p[DIM_NUM*4];
  double pn[DIM_NUM] = { 2.0, -4.0, -6.0 };
  double pp[DIM_NUM] = { -6.0, 0.0, 0.0 };
  double t[DIM_NUM*3];
  int test;
  int test_num = 6;

  cout << "\n";
  cout << "TEST030\n";
  cout << "  HALFSPACE_NORM_TRIANGLE_INT_3D finds\n";
  cout << "    intersection points of a normal form\n";
  cout << "    halfspace and a triangle.\n";

  r8vec_print ( DIM_NUM, pp, "  A point on the plane:" );
  r8vec_print ( DIM_NUM, pn, "  The normal vector:" );

  for ( test = 0; test < test_num; test++ )
  {
    if ( test == 0 )
    {
      t[0+0*3] =  0.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] = -1.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] = -2.0;
    }
    else if ( test == 1 )
    {
      t[0+0*3] = -6.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] = -1.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] = -2.0;
    }
    else if ( test == 2 )
    {
      t[0+0*3] =  0.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] =  3.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] =  2.0;
    }
    else if ( test == 3 )
    {
      t[0+0*3] = -6.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] =  4.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] =  3.0;
    }
    else if ( test == 4 )
    {
      t[0+0*3] = -8.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] = -1.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] = -2.0;
    }
    else if ( test == 5 )
    {
      t[0+0*3] =  0.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;

      t[0+1*3] =  0.0;
      t[1+1*3] =  4.0;
      t[2+1*3] =  0.0;

      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] =  4.0;
    }


    cout << "\n";
    cout << "  Case " << test << "\n";
    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    int_num = halfspace_norm_triangle_int_3d ( pp, pn, t, p );

    cout << "\n";
    cout << "  Number of intersection points is " << int_num << "\n";
    cout << "\n";

    if ( 0 < int_num )
    {
      for ( j = 0; j < int_num; j++ )
      {
        cout << "  " << setw(4)  << j
             << "  " << setw(10) << p[0+j*3]
             << "  " << setw(10) << p[1+j*3]
             << "  " << setw(10) << p[2+j*3] << "\n";
      }
    }

  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test031 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST031 tests HAVERSINE.
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
{
  double d;
  double hx;
  double pi = 3.141592653589793;
  int test;
  int test_num = 12;
  double x;

  cout << "\n";
  cout << "TEST031\n";
  cout << "  HAVERSINE computes the haversine of an angle.\n";
  cout << "\n";
  cout << "  Degrees  Radians  Haversine\n";
  cout << "\n";

  for ( test = 0; test <= test_num; test++ )
  {
    x = ( double ) ( test ) * 2.0 * pi / ( double ) ( test_num );
    d = radians_to_degrees ( x );
    hx = haversine ( x );
    cout << "  " << setw(8)  << d
         << "  " << setw(8)  << x
         << "  " << setw(12) << hx << "\n";
  }

  return;
}
//****************************************************************************80

void test0315 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0315 tests HEXAGON_CONTAINS_POINT_2D.
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
{
# define N 40
# define DIM_NUM 2

  double h[DIM_NUM*6] = {
    0.2, 0.4,
    0.4, 0.2,
    0.8, 0.0,
    1.0, 0.6,
    0.4, 1.0,
    0.2, 0.8 };
  int i;
  int j;
  double p[DIM_NUM];
  double xhi = 1.0;
  double xlo = 0.0;
  double yhi = 1.0;
  double ylo = 0.0;

  cout << "\n";
  cout << "TEST0315\n";
  cout << "  HEXAGON_CONTAINS_POINT_2D reports if a hexagon\n";
  cout << "  contains a point.\n";
  cout << "\n";
  cout << "  We will call the function repeatedly, and draw\n";
  cout << "  a sketch of an irregular hexagon in the unit square.\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    p[1] = ( ( double ) ( N - i     ) * yhi
           + ( double ) (     i - 1 ) * ylo )
           / ( double ) ( N     - 1 );

    cout << "  ";
    for ( j = 1; j <= N; j++ )
    {
      p[0] = ( ( double ) ( N - j     ) * xlo
             + ( double ) (     j - 1 ) * xhi )
             / ( double ) ( N     - 1 );

      if ( hexagon_contains_point_2d ( h, p ) )
      {
        cout << '*';
      }
      else
      {
        cout << '-';
      }
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test032 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST032 tests HEXAGON_SHAPE_2D.
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
{
# define DIM_NUM 2

  double angle;
  int i;
  double p[DIM_NUM];

  cout << "\n";
  cout << "TEST032\n";
  cout << "  HEXAGON_SHAPE_2D: points on a unit hexagon.\n";
  cout << "\n";
  cout << "  Angle    X    Y\n";
  cout << "\n";

  for ( i = -10; i <= 370; i = i + 10 )
  {
    angle = ( double ) ( i );
    hexagon_shape_2d ( angle, p );
    cout << "  " << setw(12) << angle
         << "  " << setw(12) << p[0]
         << "  " << setw(12) << p[1] << "\n";
  }

  return;
}
//****************************************************************************80

void test0321 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0321 tests HEXAGON_VERTICES_2D.
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
{
# define DIM_NUM 2

  double p[DIM_NUM*6];

  cout << "\n";
  cout << "TEST0321\n";
  cout << "  HEXAGON_VERTICES_2D: the vertices of the unit hexagon.\n";

  hexagon_vertices_2d ( p );

  r8mat_transpose_print ( DIM_NUM, 6, p, "  Vertices:" );

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0322 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0322 tests I4COL_FIND_ITEM.
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
{
# define M 5
# define N 4
# define TEST_NUM 3

  int a[M*N];
  int col;
  int i;
  int item;
  int item_test[TEST_NUM] = { 34, 12, 90 };
  int j;
  int row;
  int test;

  cout << "\n";
  cout << "TEST0322\n";
  cout << "  I4COL_FIND_ITEM finds the first occurrence of\n";
  cout << "  an item in an integer array of columns.\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = 10 * ( i + 1 ) + ( j + 1 );
    }
  }

  i4mat_print ( M, N, a, "  The matrix of columns:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    item = item_test[test];

    i4col_find_item ( M, N, a, item, &row, &col );

    cout << "  Item " << item
         << " occurs in row " << row
         << " and column " << col << "\n";
  }

  return;
# undef M
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test0323 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0323 tests I4COL_FIND_PAIR_WRAP.
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
{
# define M 5
# define N 4
# define TEST_NUM 5

  int a[M*N];
  int col;
  int i;
  int item1;
  int item1_test[TEST_NUM] = { 22, 32, 22, 54, 54 };
  int item2;
  int item2_test[TEST_NUM] = { 32, 22, 23, 14, 11 };
  int j;
  int row;
  int test;

  cout << "\n";
  cout << "TEST0323\n";
  cout << "  I4COL_FIND_PAIR_WRAP finds the first occurrence of\n";
  cout << "  a pair of item in an integer array of columns.\n";
  cout << "  Items in the array are ordered by column, and\n";
  cout << "  wraparound is allowed.\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = 10 * ( i + 1 ) + ( j + 1 );
    }
  }

  i4mat_print ( M, N, a, "  The matrix of columns:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    item1 = item1_test[test];
    item2 = item2_test[test];

    i4col_find_pair_wrap ( M, N, a, item1, item2, &row, &col );

    cout << "  Item " << item1
         << " followed by item " << item2
         << " occurs in row " << row
         << " and column " << col << "\n";
  }

  return;
# undef M
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test0325 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0325 tests ICOS_SIZE, ICOS_SHAPE, SHAPE_PRINT_3D.
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
{
# define DIM_NUM 3

  int edge_num;
  int *edge_point;
  int face_num;
  int *face_order;
  int face_order_max;
  int *face_point;
  int point_num;
  double *point_coord;

  cout << "\n";
  cout << "TEST0325\n";
  cout << "  For the icosahedron,\n";
  cout << "  ICOS_SIZE returns dimension information;\n";
  cout << "  ICOS_SHAPE returns face and order information.\n";
  cout << "  SHAPE_PRINT_3D prints this information.\n";
//
//  Get the data sizes.
//
  icos_size ( &point_num, &edge_num, &face_num, &face_order_max );

  cout << "\n";
  cout << "    Number of vertices: " << point_num << "\n";
  cout << "    Number of edges   : " << edge_num << "\n";
  cout << "    Number of faces   : " << face_num << "\n";
  cout << "    Maximum face order: " << face_order_max << "\n";
//
//  Make room for the data.
//
  point_coord = new double[DIM_NUM*point_num];
  edge_point = new int[2*edge_num];
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];
//
//  Get the data.
//
  icos_shape ( point_num, edge_num, face_num, face_order_max, point_coord,
    edge_point, face_order, face_point );
//
//  Print the data.
//
  shape_print_3d ( point_num, face_num, face_order_max,
    point_coord, face_order, face_point );

  delete [] edge_point;
  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0327 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0327 tests LINE_EXP_NORMAL_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double *normal;
  double p1[DIM_NUM] = { 1.0, 3.0 };
  double p2[DIM_NUM] = { 4.0, 0.0 };

  cout << "\n";
  cout << "TEST0327\n";
  cout << "  LINE_EXP_NORMAL_2D determines a unit normal vector\n";
  cout << "  to a given explicit line.\n";

  r8vec_print ( DIM_NUM, p1, "  Point 1: " );
  r8vec_print ( DIM_NUM, p2, "  Point 2: " );

  normal = line_exp_normal_2d ( p1, p2 );

  r8vec_print ( DIM_NUM, normal, "  Normal vector N:" );

  delete [] normal;

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test033 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST033 tests LINE_EXP_PERP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2009
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  bool flag;
  double p1[DIM_NUM] = { 1.0, 3.0 };
  double p2[DIM_NUM] = { 4.0, 0.0 };
  double *p3;
  double p3test[DIM_NUM*TEST_NUM] = {
    0.0,  0.0,
    5.0, -1.0,
    5.0,  3.0 };
  double *p4;
  int test;

  cout << "\n";
  cout << "TEST033\n";
  cout << "  LINE_EXP_PERP_2D is given an explicit line (P1,P2),\n";
  cout << "  and another point P3.  It then finds a point\n";
  cout << "  P4 on (P1,P2) so that (P1,P2) is perpendicular\n";
  cout << "  to (P3,P4).\n";

  r8vec_print ( DIM_NUM, p1, "  Point P1:" );
  r8vec_print ( DIM_NUM, p2, "  Point P2:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p3 = p3test + test * DIM_NUM;

    r8vec_print ( DIM_NUM, p3, "  Point P3:" );

    p4 = line_exp_perp_2d ( p1, p2, p3, &flag );

    r8vec_print ( DIM_NUM, p4, "  Point P4:" );

    delete [] p4;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0335 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0335 tests LINE_EXP_POINT_DIST_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double dist;
  double *p;
  double p1[DIM_NUM] = { 1.0, 3.0 };
  double p2[DIM_NUM] = { 4.0, 0.0 };
  double p_test[DIM_NUM*TEST_NUM] = {
    0.0,  0.0,
    5.0, -1.0,
    5.0,  3.0 };
  int test;

  cout << "\n";
  cout << "TEST0335\n";
  cout << "  LINE_EXP_POINT_DIST_2D finds the distance from\n";
  cout << "    an explicit line to a point in 2D.\n";

  r8vec_print ( DIM_NUM, p1, "  Point 1: " );
  r8vec_print ( DIM_NUM, p2, "  Point 2: " );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    r8vec_print ( DIM_NUM, p, "  Point: " );

    dist = line_exp_point_dist_2d ( p1, p2, p );

    cout << "  Distance = " << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0336 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0336 tests LINE_EXP_POINT_DIST_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 3

  double dist;
  double *p;
  double p1[DIM_NUM] = { 1.0, 3.0, 2.0 };
  double p2[DIM_NUM] = { 4.0, 0.0, 1.0 };
  double p_test[DIM_NUM*TEST_NUM] = {
    0.0,  0.0, 2.0,
    5.0, -1.0, 1.0,
    5.0,  3.0, 3.0 };
  int test;

  cout << "\n";
  cout << "TEST0336\n";
  cout << "  LINE_EXP_POINT_DIST_3D finds the distance\n";
  cout << "    from an explicit line to a point in 3D.\n";

  r8vec_print ( DIM_NUM, p1, "  Point 1:" );
  r8vec_print ( DIM_NUM, p2, "  Point 2:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    r8vec_print ( DIM_NUM, p, "  Point:" );

    dist = line_exp_point_dist_3d ( p1, p2, p );

    cout << "  Distance = " << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0337 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0337 tests LINE_EXP_POINT_DIST_SIGNED_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double dist;
  double *p;
  double p1[DIM_NUM] = { 1.0, 3.0 };
  double p2[DIM_NUM] = { 4.0, 0.0 };
  double p_test[DIM_NUM*TEST_NUM] = {
    0.0,  0.0,
    5.0, -1.0,
    5.0,  3.0 };
  int test;

  cout << "\n";
  cout << "TEST0337\n";
  cout << "  LINE_EXP_POINT_DIST_SIGNED_2D finds the signed\n";
  cout << "    distance to a point from an explicit line.\n";

  r8vec_print ( DIM_NUM, p1, "  Point 1:" );
  r8vec_print ( DIM_NUM, p2, "  Point 2:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    r8vec_print ( DIM_NUM, p, "  Point:" );

    dist = line_exp_point_dist_signed_2d ( p1, p2, p );

    cout << "  Signed distance = " << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test034 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST034 tests LINE_EXP_POINT_NEAR_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double dist;
  double *p;
  double p1[DIM_NUM] = { 1.0, 3.0 };
  double p2[DIM_NUM] = { 4.0, 0.0 };
  double pn[DIM_NUM];
  double p_test[DIM_NUM*TEST_NUM] = {
    0.0,  0.0,
    5.0, -1.0,
    5.0,  3.0 };
  double t;
  int test;

  cout << "\n";
  cout << "TEST034\n";
  cout << "  LINE_EXP_POINT_NEAR_2D finds the point on\n";
  cout << "    a line nearest in point in 2D.\n";

  r8vec_print ( DIM_NUM, p1, "  The point P1:" );
  r8vec_print ( DIM_NUM, p2, "  The point P2:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    r8vec_print ( DIM_NUM, p, "  The point P:" );

    line_exp_point_near_2d ( p1, p2, p, pn, &dist, &t );

    r8vec_print ( DIM_NUM, pn, "  Nearest point PN:" );

    cout << "  Distance = " << dist << "\n";
    cout << "  Relative line position T = " << t << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0345 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0345 tests LINE_EXP2IMP_2D and LINE_IMP2EXP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double a1 = 1.0;
  double a2;
  double b1 = 2.0;
  double b2;
  double c1 = 3.0;
  double c2;
  double p1[DIM_NUM];
  double p2[DIM_NUM];

  cout << "\n";
  cout << "TEST0345\n";
  cout << "  LINE_EXP2IMP_2D converts explicit to implicit lines.\n";
  cout << "  LINE_IMP2EXP_2D converts implicit to explicit lines.\n";

  cout << "\n";
  cout << "  Implicit line A = " << a1
       << " B = " << b1
       << " C = " << c1 << "\n";

  line_imp2exp_2d ( a1, b1, c1, p1, p2 );

  r8vec_print ( DIM_NUM, p1, "  The point P1:" );
  r8vec_print ( DIM_NUM, p2, "  The point P2:" );

  line_exp2imp_2d ( p1, p2, &a2, &b2, &c2 );

  cout << "\n";
  cout << "  Recovered implicit line A = " << a2
       << " B = " << b2
       << " C = " << c2 << "\n";

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0346 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0346 tests LINE_EXP2PAR_2D and LINE_PAR2EXP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double f1 = 1.0;
  double f2;
  double g1 = 2.0;
  double g2;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double x1 = 3.0;
  double x2;
  double y1 = 4.0;
  double y2;

  cout << "\n";
  cout << "TEST0346\n";
  cout << "  LINE_EXP2PAR_2D converts explicit to parametric lines.\n";
  cout << "  LINE_PAR2EXP_2D converts parametric to explicit lines.\n";

  cout << "\n";
  cout << "  Parametric line:\n";
  cout << "    F  = " << f1 << " G  = " << g1 << "\n";
  cout << "    X0 = " << x1 << " Y0 = " << y1 << "\n";

  line_par2exp_2d ( f1, g1, x1, y1, p1, p2 );

  r8vec_print ( DIM_NUM, p1, "  The point P1:" );
  r8vec_print ( DIM_NUM, p2, "  The point P2:" );

  line_exp2par_2d ( p1, p2, &f2, &g2, &x2, &y2 );

  cout << "\n";
  cout << "  Recovered parametric line:\n";
  cout << "    F  = " << f2 << " G  = " << g2 << "\n";
  cout << "    X0 = " << x2 << " Y0 = " << y2 << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test035 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST035 tests LINE_IMP_POINT_DIST_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double a;
  double atest[TEST_NUM] = { 2.0, 2.0, 2.0 };
  double b;
  double btest[TEST_NUM] = { 5.0, 5.0, 5.0 };
  double c;
  double ctest[TEST_NUM] = { 3.0, 3.0, 3.0 };
  double dist;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
    0.0, 6.0,
    0.0, 5.0,
    0.0, 4.0 };
  int test;

  cout << "\n";
  cout << "TEST035\n";
  cout << "  LINE_IMP_POINT_DIST_2D finds the distance from\n";
  cout << "    a point P to a line A * X + B * Y + C = 0.\n";
  cout << "\n";
  cout << "   X       Y       A       B       C       DIST\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = atest[test];
    b = btest[test];
    c = ctest[test];
    p = p_test + test * DIM_NUM;

    dist = line_imp_point_dist_2d ( a, b, c, p );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  " << setw(8) << a
         << "  " << setw(8) << b
         << "  " << setw(8) << c
         << "  " << setw(8) << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0351 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0351 tests LINE_PAR_POINT_NEAR_2D and LINE_PAR_POINT_DIST_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2013
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double dist;
  double f;
  double g;
  int i;
  double p[2];
  double *pn;
  double p_test[2*TEST_NUM] = {
    0.0,  0.0,
    5.0, -1.0,
    5.0,  3.0 };
  double t;
  int test;
  int test_num = TEST_NUM;
  double x0;
  double y0;

  cout << "\n";
  cout << "TEST0351\n";
  cout << "  LINE_PAR_POINT_NEAR_2D finds the point on\n";
  cout << "  a parametric line (X0,Y0,F,G) nearest a point P in 2D.\n";

  x0 = 1.0;
  y0 = 3.0;
  f = +1.0;
  g = -1.0;

  cout << "\n";
  cout << "  Parametric line:\n";
  cout << "  X(t) = " << x0 << " + " << f << " * t\n";
  cout << "  Y(t) = " << y0 << " + " << g << " * t\n";

  for ( test = 0; test < test_num; test++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      p[i] = p_test[i+test*2];
    }

    r8vec_print ( 2, p, "  The point P:" );

    dist = line_par_point_dist_2d ( f, g, x0, y0, p );

    cout << "  Distance = " << dist << "\n";

    pn = line_par_point_near_2d ( f, g, x0, y0, p );

    r8vec_print ( 2, pn, "  Nearest point PN:" );

    dist = r8vec_norm_affine ( 2, p, pn );

    cout << "  Distance recomputed = " << dist << "\n";

    delete [] pn;
  }
  return;
# undef TEST_NUM
}
//****************************************************************************80

void test0352 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0352 tests LINE_PAR_POINT_DIST_3D and LINE_PAR_POINT_NEAR_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 April 2013
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double dist;
  double f;
  double g;
  double h;
  int i;
  double p[3];
  double *pn;
  double p_test[3*TEST_NUM] = {
    0.0,  0.0, 2.0, 
    5.0, -1.0, 1.0, 
    5.0,  3.0, 3.0 };
  int test;
  int test_num = TEST_NUM;
  double x0;
  double y0;
  double z0;

  cout << "\n";
  cout << "TEST0352\n";
  cout << "  LINE_PAR_POINT_DIST_3D finds the distance\n";
  cout << "  from a parametric line to a point in 3D.\n";

  x0 = 1.0;
  y0 = 3.0;
  z0 = 2.0;

  f = +3.0;
  g = -3.0;
  h = -1.0;

  cout << "\n";
  cout << "  Parametric line:\n";
  cout << "  X(t) = " << x0 << " + " << f << " * t\n";
  cout << "  Y(t) = " << y0 << " + " << g << " * t\n";
  cout << "  Z(t) = " << z0 << " + " << h << " * t\n";

  for ( test = 0; test < test_num; test++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      p[i] = p_test[i+test*3];
    }

    r8vec_print ( 3, p, "  The point P:" );

    dist = line_par_point_dist_3d ( f, g, h, x0, y0, z0, p );

    cout << "  Distance = " << dist << "\n";

    pn = line_par_point_near_3d ( f, g, h, x0, y0, z0, p );

    r8vec_print ( 3, pn, "  Nearest point PN:" );

    dist = r8vec_norm_affine ( 3, p, pn );

    cout << "  Distance recomputed = " << dist << "\n";

    delete [] pn;
  }
  return;
# undef TEST_NUM
}
//****************************************************************************80

void test038 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST038 tests LINES_EXP_ANGLE_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  double angle;
  double *p1;
  double p1_test[DIM_NUM*TEST_NUM] = {
    0.0, 0.0, 0.0,
    1.0, 2.0, 0.0 };
  double *p2;
  double p2_test[DIM_NUM*TEST_NUM] = {
    1.0, 2.0, 0.0,
    1.0, 2.0, 0.0 };
  double *q1;
  double q1_test[DIM_NUM*TEST_NUM] = {
    0.0, 3.0,  3.0,
    1.0, 2.0, -1.0 };
  double *q2;
  double q2_test[DIM_NUM*TEST_NUM] = {
    3.0, 0.0, 3.0,
    1.0, 2.0, 3.0 };
  int test;

  cout << "\n";
  cout << "TEST038\n";
  cout << "  LINES_EXP_ANGLE_3D finds the angle between\n";
  cout << "    two explicit lines in 3D;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p1_test + test * DIM_NUM;
    p2 = p2_test + test * DIM_NUM;
    q1 = q1_test + test * DIM_NUM;
    q2 = q2_test + test * DIM_NUM;

    angle = lines_exp_angle_3d ( p1, p2, q1, q2 );

    cout << "\n";
    cout << "  Angle between lines is " << angle << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80*

void test0385 ( )

//****************************************************************************80*
//
//  Purpose:
//
//    TEST0385 tests LINES_EXP_DIST_3D and LINES_EXP_DIST_3D_2;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  double dist;
  double dist2;
  int i;
  double *p1;
  double p1_test[DIM_NUM*TEST_NUM] = {
    0.0,  0.0, 0.0,
    4.0, -3.0, 0.0 };
  double *p2;
  double p2_test[DIM_NUM*TEST_NUM] = {
    1.0, 2.0, 0.0,
   -8.0, 6.0, 0.0 };
  double *q1;
  double q1_test[DIM_NUM*TEST_NUM] = {
    0.0, 3.0,  3.0,
    3.0, 4.0, -1.0 };
  double *q2;
  double q2_test[DIM_NUM*TEST_NUM] = {
    3.0, 0.0, 3.0,
    3.0, 4.0, 3.0 };
  int test;

  cout << "\n";
  cout << "TEST0385\n";
  cout << "  LINES_EXP_DIST_3D finds the distance between\n";
  cout << "    two explicit lines in 3D;\n";
  cout << "  LINES_EXP_DIST_3D_2 finds the distance between\n";
  cout << "    two explicit lines in 3D;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p1_test + test * DIM_NUM;
    p2 = p2_test + test * DIM_NUM;
    q1 = q1_test + test * DIM_NUM;
    q2 = q2_test + test * DIM_NUM;

    dist = lines_exp_dist_3d ( p1, p2, q1, q2 );
    dist2 = lines_exp_dist_3d_2 ( p1, p2, q1, q2 );

    cout << "\n";
    cout << "\n";
    cout << "  P1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p1[i];
    }
    cout << "\n";
    cout << "  P2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p2[i];
    }
    cout << "\n";
    cout << "  Q1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q1[i];
    }
    cout << "\n";
    cout << "  Q2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q2[i];
    }
    cout << "\n";
    cout << "  LINES_EXP_DIST_3D =   " << dist << "\n";
    cout << "  LINES_EXP_DIST_3D_2 = " << dist2 << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80*

void test03855 ( )

//****************************************************************************80*
//
//  Purpose:
//
//    TEST03855 tests LINES_EXP_DIST_3D and LINES_EXP_DIST_3D_2;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  int i;
  double *p1;
  double p1_test[DIM_NUM*TEST_NUM] = {
    0.0,  0.0, 0.0,
    4.0, -3.0, 0.0 };
  double *p2;
  double p2_test[DIM_NUM*TEST_NUM] = {
    1.0, 2.0, 0.0,
   -8.0, 6.0, 0.0 };
  double pn[DIM_NUM];
  double *q1;
  double q1_test[DIM_NUM*TEST_NUM] = {
    0.0, 3.0,  3.0,
    3.0, 4.0, -1.0 };
  double *q2;
  double q2_test[DIM_NUM*TEST_NUM] = {
    3.0, 0.0, 3.0,
    3.0, 4.0, 3.0 };
  double qn[DIM_NUM];
  int test;

  cout << "\n";
  cout << "TEST03855\n";
  cout << "  LINES_EXP_NEAR_3D finds nearest points on\n";
  cout << "    two explicit lines in 3D;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p1_test + test * DIM_NUM;
    p2 = p2_test + test * DIM_NUM;
    q1 = q1_test + test * DIM_NUM;
    q2 = q2_test + test * DIM_NUM;

    cout << "\n";
    cout << "\n";
    cout << "  P1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p1[i];
    }
    cout << "\n";
    cout << "  P2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p2[i];
    }
    cout << "\n";
    cout << "  Q1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q1[i];
    }
    cout << "\n";
    cout << "  Q2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q2[i];
    }
    cout << "\n";

    lines_exp_near_3d ( p1, p2, q1, q2, pn, qn );

    cout << "\n";
    cout << "  PN:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << pn[i];
    }
    cout << "\n";
    cout << "  QN:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << qn[i];
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80*

void test0386 ( )

//****************************************************************************80*
//
//  Purpose:
//
//    TEST0386 tests LINES_EXP_EQUAL_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 6

  bool equal;
  int i;
  double *p1;
  double p1_test[DIM_NUM*TEST_NUM] = {
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0 };
  double *p2;
  double p2_test[DIM_NUM*TEST_NUM] = {
    1.0, 2.0,
    1.0, 2.0,
    1.0, 2.0,
    1.0, 2.0,
    1.0, 2.0,
    1.0, 2.0 };
  double *q1;
  double q1_test[DIM_NUM*TEST_NUM] = {
    0.0,  0.0,
    1.0,  2.0,
    0.0,  0.0,
    7.0, 14.0,
    1.0,  2.0,
    0.0, 10.0 };
  double *q2;
  double q2_test[DIM_NUM*TEST_NUM] = {
    1.0,  2.0,
    0.0,  0.0,
    2.0,  4.0,
    5.5, 11.0,
    3.0,  5.0,
    1.0, 12.0 };
  int test;

  cout << "\n";
  cout << "TEST0386\n";
  cout << "  LINES_EXP_EQUAL_2D tries to determine if two\n";
  cout << "    explicit lines in 2D are equal;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p1_test + test * DIM_NUM;
    p2 = p2_test + test * DIM_NUM;
    q1 = q1_test + test * DIM_NUM;
    q2 = q2_test + test * DIM_NUM;

    cout << "\n";
    cout << "  P1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p1[i];
    }
    cout << "\n";
    cout << "  P2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p2[i];
    }
    cout << "\n";
    cout << "  Q1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q1[i];
    }
    cout << "\n";
    cout << "  Q2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q2[i];
    }
    cout << "\n";

    equal = lines_exp_equal_2d ( p1, p2, q1, q2 );

    if ( equal )
    {
      cout << "  The lines are equal.\n";
    }
    else
    {
      cout << "  The lines are distinct.\n";
    }
  }
  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80*

void test039 ( )

//****************************************************************************80*
//
//  Purpose:
//
//    TEST039 tests LINES_EXP_INT_2D;
//
//  Discussion:
//
//    Test #1:
//
//      x + 2y -  4 = 0
//      x -  y -  1 = 0
//
//    Test #2:
//
//      x + 2y -  4 = 0
//     2x + 4y -  1 = 0
//
//    Test #3:
//
//      x + 2y -  4 = 0
//    -3x - 6y + 12 = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  int i;
  int ival;
  double p[DIM_NUM];
  double *p1;
  double p1_test[DIM_NUM*TEST_NUM] = {
    0.0, 0.0,
    0.0, 2.0,
    0.0, 2.0 };
  double *p2;
  double p2_test[DIM_NUM*TEST_NUM] = {
    4.0, 0.0,
    4.0, 0.0,
    4.0, 0.0 };
  double *q1;
  double q1_test[DIM_NUM*TEST_NUM] = {
    0.0, -1.0,
    0.0,  0.25,
    0.0,  2.0 };
  double *q2;
  double q2_test[DIM_NUM*TEST_NUM] = {
    1.0,  0.0,
    0.5,  0.0,
    4.0,  0.0 };
  int test;

  cout << "\n";
  cout << "TEST039\n";
  cout << "  LINES_EXP_INT_2D finds intersections of\n";
  cout << "    two explicit lines in 2D;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p1_test + test * DIM_NUM;
    p2 = p2_test + test * DIM_NUM;
    q1 = q1_test + test * DIM_NUM;
    q2 = q2_test + test * DIM_NUM;

    cout << "\n";
    cout << "  P1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p1[i];
    }
    cout << "\n";
    cout << "  P2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p2[i];
    }
    cout << "\n";
    cout << "  Q1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q1[i];
    }
    cout << "\n";
    cout << "  Q2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q2[i];
    }
    cout << "\n";

    lines_exp_int_2d ( p1, p2, q1, q2, &ival, p );

    if ( ival == 1 )
    {
      cout << "  Intersection at";
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << "  " << setw(12) << p[i];
      }
      cout << "\n";
    }
    else if ( ival == 0 )
    {
      cout << "  Lines are parallel, no intersection.\n";
    }
    else if ( ival == 2 )
    {
      cout << "  Lines are coincident.\n";
    }
    else
    {
      cout << "  Unknown return value of IVAL = " << ival << "\n";
    }
  }
  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test040 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST040 tests LINES_IMP_ANGLE_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double a1;
  double a2;
  double angle;
  double b1;
  double b2;
  double c1;
  double c2;

  cout << "\n";
  cout << "TEST040\n";
  cout << "  For two lines written in implicit form:\n";
  cout << "\n";
  cout << "  LINES_IMP_ANGLE_2D finds the angle;\n";
  cout << "\n";
//
//  x + 2y - 4 = 0
//
  a1 =   1.0;
  b1 =   2.0;
  c1 = - 4.0;

  cout << "\n";
  cout << "  Line 1 coefficients: " << a1 << "  " << b1 << "  " << c1 << "\n";
//
//  x - y - 1 = 0
//
  a2 =   1.0;
  b2 = - 1.0;
  c2 = - 1.0;
  cout << "  Line 2 coefficients: " << a2 << "  " << b2 << "  " << c2 << "\n";

  angle = lines_imp_angle_2d ( a1, b1, c1, a2, b2, c2 );

  cout << "\n";
  cout << "  Angle between lines is " << radians_to_degrees ( angle ) << "\n";

  cout << "\n";
  cout << "  Line 1 coefficients: " << a1 << "  " << b1 << "  " << c1 << "\n";
//
//  2x + 4y - 1 = 0
//
  a2 =   2.0;
  b2 = + 4.0;
  c2 = - 1.0;
  cout << "  Line 2 coefficients: " << a2 << "  " << b2 << "  " << c2 << "\n";

  angle = lines_imp_angle_2d ( a1, b1, c1, a2, b2, c2 );

  cout << "\n";
  cout << "  Angle between lines is " << radians_to_degrees ( angle ) << "\n";

  cout << "\n";
  cout << "  Line 1 coefficients: " << a1 << "  " << b1 << "  " << c1 << "\n";
//
//  -3x - 6y +12 = 0
//
  a2 =  - 3.0;
  b2 =  - 6.0;
  c2 = + 12.0;
  cout << "  Line 2 coefficients: " << a2 << "  " << b2 << "  " << c2 << "\n";

  angle = lines_imp_angle_2d ( a1, b1, c1, a2, b2, c2 );

  cout << "\n";
  cout << "  Angle between lines is " << radians_to_degrees ( angle ) << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test041 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST041 tests LINES_IMP_DIST_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double a1;
  double a1_test[TEST_NUM] = {  4.0,  2.0, 1.0 };
  double a2;
  double a2_test[TEST_NUM] = {  4.0,  4.0, 2.0 };
  double b1;
  double b1_test[TEST_NUM] = { -1.0, -1.0, 2.0 };
  double b2;
  double b2_test[TEST_NUM] = { -1.0, -2.0, 3.0 };
  double c1;
  double c1_test[TEST_NUM] = {  3.0,  0.0, 2.0 };
  double c2;
  double c2_test[TEST_NUM] = { 12.0,  6.0, 1.0 };
  double dist;
  int test;

  cout << "\n";
  cout << "TEST041\n";
  cout << "  LINES_IMP_DIST_3D finds the distance between\n";
  cout << "    two implicit lines in 2D.\n";
  cout << "\n";
  cout << "   A1      B1      C1      A2      B2      C2   DIST\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a1 = a1_test[test];
    b1 = b1_test[test];
    c1 = c1_test[test];
    a2 = a2_test[test];
    b2 = b2_test[test];
    c2 = c2_test[test];

    dist = lines_imp_dist_2d ( a1, b1, c1, a2, b2, c2 );

    cout << "  " << setw(8) << a1
         << "  " << setw(8) << b1
         << "  " << setw(8) << c1
         << "  " << setw(8) << a2
         << "  " << setw(8) << b2
         << "  " << setw(8) << c2
         << "  " << setw(8) << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0415 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0415 tests LINES_IMP_INT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double a1;
  double a2;
  double b1;
  double b2;
  double c1;
  double c2;
  int ival;
  double p[2];

  cout << "\n";
  cout << "TEST0415\n";
  cout << "  For two lines written in implicit form:\n";
  cout << "\n";
  cout << "  LINES_IMP_INT_2D finds the intersection.\n";
  cout << "\n";
//
//  x + 2y - 4 = 0
//
  a1 =   1.0;
  b1 =   2.0;
  c1 = - 4.0;

  cout << "\n";
  cout << "  Line 1 coefficients: " << a1 << "  " << b1 << "  " << c1 << "\n";
//
//  x - y - 1 = 0
//
  a2 =   1.0;
  b2 = - 1.0;
  c2 = - 1.0;
  cout << "  Line 2 coefficients: " << a2 << "  " << b2 << "  " << c2 << "\n";

  lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, &ival, p );

  if ( ival == 1 )
  {
    cout << "  Intersection at " << p[0] << "  " << p[1] << "\n";
  }
  else if ( ival == 0 )
  {
    cout << "  Lines are parallel, no intersection.\n";
  }
  else if ( ival == 2 )
  {
    cout << "  Lines are coincident.\n";
  }
  else
  {
    cout << "  Unknown return value of ival = " << ival << "\n";
  }

  cout << "\n";
  cout << "  Line 1 coefficients: " << a1 << "  " << b1 << "  " << c1 << "\n";
//
//  2x + 4y - 1 = 0
//
  a2 =   2.0;
  b2 = + 4.0;
  c2 = - 1.0;
  cout << "  Line 2 coefficients: " << a2 << "  " << b2 << "  " << c2 << "\n";

  lines_imp_int_2d (a1, b1, c1, a2, b2, c2, &ival, p );

  if ( ival == 1 )
  {
    cout << "  Intersection at " << p[0] << "  " << p[1] << "\n";
  }
  else if ( ival == 0 )
  {
    cout << "  Lines are parallel, no intersection.\n";
  }
  else if ( ival == 2 )
  {
    cout << "  Lines are coincident.\n";
  }
  else
  {
    cout << "  Unknown return value of IVAL = " << ival << "\n";
  }

  cout << "\n";
  cout << "  Line 1 coefficients: " << a1 << "  " << b1 << "  " << c1 << "\n";
//
//  -3x - 6y +12 = 0
//
  a2 =  - 3.0;
  b2 =  - 6.0;
  c2 = + 12.0;
  cout << "  Line 2 coefficients: " << a2 << "  " << b2 << "  " << c2 << "\n";

  lines_imp_int_2d  ( a1, b1, c1, a2, b2, c2, &ival, p );

  if ( ival == 1 )
  {
    cout << "  Intersection at " << p[0] << "  " << p[1] << "\n";
  }
  else if ( ival == 0 )
  {
    cout << "  Lines are parallel, no intersection.\n";
  }
  else if ( ival == 2 )
  {
    cout << "  Lines are coincident.\n";
  }
  else
  {
    cout << "  Unknown return value of IVAL = " << ival << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0416 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0416 tests LINES_PAR_INT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double f1;
  double f2;
  double g1;
  double g2;
  double pint[DIM_NUM];
  double t1;
  double t2;
  double x1;
  double x2;
  double y1;
  double y2;

  cout << "\n";
  cout << "TEST0416\n";
  cout << "  LINES_PAR_INT_2D finds the intersection of\n";
  cout << "  two lines written in parametric form.\n";
  cout << "\n";
//
//  x - 2y = -1
//
  x1 =  0.0;
  y1 =  1.0;
  f1 =  2.0;
  g1 =  1.0;

  cout << "\n";
  cout << "  Line 1 parameters:" << setw(8) << x1
       << "  " << setw(8) << y1
       << "  " << setw(8) << f1
       << "  " << setw(8) << g1 << "\n";
//
//  x + y - 8 = 0
//
  x2 = 10.0;
  y2 = -2.0;
  f2 =  1.0;
  g2 =  1.0;

  cout << "\n";
  cout << "  Line 2 parameters:" << setw(8) << x2
       << "  " << setw(8) << y2
       << "  " << setw(8) << f2
       << "  " << setw(8) << g2 << "\n";

  lines_par_int_2d ( f1, g1, x1, y1, f2, g2, x2, y2, &t1, &t2, pint );

  cout << "\n";
  cout << "  Line 1 evaluated at T1:\n";
  cout << "\n";
  cout << "    T1 =   " << t1 << "\n";
  cout << "    X(T1)= " << x1 + f1 * t1 << "\n";
  cout << "    Y(T1)= " << y1 + g1 * t1 << "\n";
  cout << "\n";
  cout << "  Line 2 evaluated at T2:\n";
  cout << "\n";
  cout << "    T2 =   " << t2 << "\n";
  cout << "    X(T2)= " << x2 + f2 * t2 << "\n";
  cout << "    Y(T2)= " << y2 + g2 * t2 << "\n";
  cout << "\n";
  cout << "  Reported intersection PINT:\n";
  cout << "\n";
  cout << "    " << pint[0] << "\n";
  cout << "    " << pint[1] << "\n";

  return;
}
//****************************************************************************80

void test0418 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0418 tests SEGMENTS_CURVATURE_2D.
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
{
# define DIM_NUM 2
# define TEST_NUM 13

  double curvature;
  double p1[DIM_NUM] = { 0.0, 0.0 };
  double p2[DIM_NUM] = { 1.0, 0.0 };
  double p3[DIM_NUM];
  double pi = 3.141592653589793;
  double theta;
  double theta_degrees;
  int test;

  cout << "\n";
  cout << "TEST0418\n";
  cout << "  SEGMENTS_CURVATURE_2D computes the local curvature\n";
  cout << "  defined by the line segments [P1,P2] and [P2,P3].\n";

  cout << "\n";
  cout << "  Our three points are:\n";
  cout << "\n";
  cout << "    P1 = (0,0)\n";
  cout << "    P2 = (1,0)\n";
  cout << "    P3 = (C,S)\n";
  cout << "\n";
  cout << "  C = cosine ( theta), S = sine ( theta ).\n";
  cout << "\n";
  cout << "  Test  Theta  Curvature\n";
  cout << "\n";

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    theta = 2.0 * pi * ( double ) ( test - 1 )
      / ( double ) ( TEST_NUM - 1 );

    theta_degrees = 360.0 * ( double ) ( test - 1 )
      / ( double ) ( TEST_NUM - 1 );

    p3[0] = cos ( theta );
    p3[1] = sin ( theta );

    curvature = segments_curvature_2d ( p1, p2, p3 );

    cout << "  " << setw(4)  << test
         << "  " << setw(5)  << theta_degrees
         << "  " << setw(14) << curvature << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test042 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST042 tests SEGMENTS_DIST_2D.
//
//  Discussion:
//
//    Case 1, parallel, not coincident.
//    Case 2, parallel, coincident, overlapping.
//    Case 3, parallel, coincident, disjoint.
//    Case 4, nonparallel, intersecting.
//    Case 5, nonparallel, disjoint.
//    Case 6 and 7, should be same, because simply a translation by 50;
//    Case 8 and 9, answers should be same.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 9

  double dist;
  int i;
  double *p1;
  double p1_test[DIM_NUM*TEST_NUM] = {
    2.0,   3.0,
    2.0,   3.0,
    2.0,   3.0,
    2.0,   3.0,
    2.0,   3.0,
   57.0,  53.0,
    7.0,   3.0,
    0.0,   0.0,
  -10.0, -10.0 };
  double *p2;
  double p2_test[DIM_NUM*TEST_NUM] = {
    8.0,   6.0,
    8.0,   6.0,
    8.0,   6.0,
    8.0,   6.0,
    8.0,   6.0,
   58.0,  53.0,
    8.0,   3.0,
  100.0, 100.0,
  100.0, 100.0 };
  double *q1;
  double q1_test[DIM_NUM*TEST_NUM] = {
    8.0,  3.0,
    4.0,  4.0,
   14.0,  9.0,
    0.0,  8.0,
    7.0,  3.0,
   65.0, 45.0,
   15.0, -5.0,
   50.0,  0.0,
   50.0,  0.0 };
  double *q2;
  double q2_test[DIM_NUM*TEST_NUM] = {
   14.0,  6.0,
   14.0,  9.0,
   16.0, 10.0,
    5.0,  3.0,
    9.0, -1.0,
   57.0, 53.0,
    7.0,  3.0,
   60.0,  0.0,
   60.0,  0.0 };
  int test;

  cout << "\n";
  cout << "TEST042\n";
  cout << "  SEGMENTS_DIST_2D computes the distance between\n";
  cout << "    line segments in 2D.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p1_test + test * DIM_NUM;
    p2 = p2_test + test * DIM_NUM;
    q1 = q1_test + test * DIM_NUM;
    q2 = q2_test + test * DIM_NUM;

    dist = segments_dist_2d ( p1, p2, q1, q2 );

    cout << "\n";

    if ( test == 0 )
    {
      cout << "  Same slope, different intercepts.\n";
    }
    else if ( test == 1 )
    {
      cout << "  Same slope, same intercepts, overlapping.\n";
      cout << "  Distance should be 0.\n";
    }
    else if ( test == 2 )
    {
      cout << "  Same slope, same intercepts, disjoint.\n";
      cout << "  Distance should be sqrt(45)=6.7082038\n";
    }
    else if ( test == 3 )
    {
      cout << "  Different slopes, intersecting.\n";
      cout << "  Distance should be 0.\n";
    }
    else if ( test == 4 )
    {
      cout << "  Different slopes, not intersecting.\n";
    }
    else if ( test == 5 )
    {
      cout << "  Simple problem.\n";
      cout << "  Distance should be 0.\n";
    }
    else if ( test == 6 )
    {
      cout << "  Same data, translated by 50.\n";
      cout << "  Distance should be 0.\n";
    }
    else if ( test == 7 )
    {
      cout << "  Diagonal and horizontal.\n";
      cout << "  Distance should be sqrt(2500/2)=35.355339\n";
    }
    else if ( test == 8 )
    {
      cout << "  Same data, except first segment extended.\n";
      cout << "  Distance should be sqrt(2500/2)=35.355339\n";
    }
    cout << "\n";
    cout << "  P1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p1[i];
    }
    cout << "\n";
    cout << "  P2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p2[i];
    }
    cout << "  Q1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q1[i];
    }
    cout << "\n";
    cout << "  Q2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q2[i];
    }
    cout << "  Distance = " << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test043 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST043 tests SEGMENTS_DIST_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double dist;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];
  double p4[DIM_NUM];

  cout << "\n";
  cout << "TEST043\n";
  cout << "  SEGMENTS_DIST_3D computes the distance between two\n";
  cout << "  line segments in 3D.\n";
  cout << " \n";
  cout << "  Case   Computed    True\n";
  cout << " \n";
//
//  Case 0, colinear, not coincident.
//
//  LS1: (2,3,0) + t * (2,1,0) for t = 0 to 3.
//  LS2: (11,6,4) + t * (2,1,0) for t = 0 to 3.
//  Distance is 5.
//
  p1[0] = 0.0;
  p1[1] = 0.0;
  p1[2] = 0.0;

  p2[0] = 1.0;
  p2[1] = 0.0;
  p2[2] = 0.0;

  p3[0] = 2.0;
  p3[1] = 0.0;
  p3[2] = 0.0;

  p4[0] = 4.0;
  p4[1] = 0.0;
  p4[2] = 0.0;

  dist = segments_dist_3d ( p1, p2, p3, p4 );

  cout << "  " << 0 << "  " << dist << "  " << 1.0 << "\n";
//
//  Case 1, parallel, not coincident.
//
//  LS1: (2,3,0) + t * (2,1,0) for t = 0 to 3.
//  LS2: (11,6,4) + t * (2,1,0) for t = 0 to 3.
//  Distance is 5.
//
  p1[0] = 2.0;
  p1[1] = 3.0;
  p1[2] = 0.0;

  p2[0] = 8.0;
  p2[1] = 6.0;
  p2[2] = 0.0;

  p3[0] = 11.0;
  p3[1] =  6.0;
  p3[2] =  4.0;

  p4[0] = 17.0;
  p4[1] =  9.0;
  p4[2] =  4.0;

  dist = segments_dist_3d ( p1, p2, p3, p4 );

  cout << "  " << setw(3) << 1
       << "  " << setw(12) << dist << "  " << 5.0 << "\n";
//
//  Case 2, parallel, coincident, overlapping.
//
//  (1,2,3) + t * ( 1,-1,2)
//  LS1: t = 0 to t = 3.
//  Distance is 0.
//
  p1[0] = 1.0;
  p1[1] = 2.0;
  p1[2] = 3.0;

  p2[0] =  4.0;
  p2[1] = -1.0;
  p2[2] =  9.0;

  p3[0] = 3.0;
  p3[1] = 0.0;
  p3[2] = 7.0;

  p4[0] =  6.0;
  p4[1] = -3.0;
  p4[2] = 13.0;

  dist = segments_dist_3d ( p1, p2, p3, p4 );

  cout << "  " << setw(3) << 2
       << "  " << setw(12) << dist << "  " << 0.0 << "\n";
//
//  Case 3, parallel, coincident, disjoint.
//
//  LS1: (3,4,5) + t * ( 2,2,1) for 0 <= t <= 2.
//  LS2: (3,4,5) + t * ( 2,2,1) for 3 <= t <= 5.
//  Distance = 3.
//
  p1[0] = 3.0;
  p1[1] = 4.0;
  p1[2] = 5.0;

  p2[0] =  7.0;
  p2[1] =  8.0;
  p2[2] =  7.0;

  p3[0] =  9.0;
  p3[1] = 10.0;
  p3[2] =  8.0;

  p4[0] =  13.0;
  p4[1] =  14.0;
  p4[2] =  10.0;

  dist = segments_dist_3d ( p1, p2, p3, p4 );

  cout << "  " << setw(3) << 3
       << "  " << setw(12) << dist << "  " << 3.0 << "\n";
//
//  Case 4, nonparallel, could intersect, and does intersect.
//
//  L1: (1,1,1) + t * (0,1,2)
//  L2: (0,2,3) + t * (1,0,0)
//  intersect at (1,2,3)
//  Distance is 0.
//
  p1[0] = 1.0;
  p1[1] = 1.0;
  p1[2] = 1.0;

  p2[0] = 1.0;
  p2[1] = 4.0;
  p2[2] = 7.0;

  p3[0] =  0.0;
  p3[1] =  2.0;
  p3[2] =  3.0;

  p4[0] =  5.0;
  p4[1] =  2.0;
  p4[2] =  3.0;

  dist = segments_dist_3d ( p1, p2, p3, p4 );

  cout << "  " << setw(3) << 4
       << "  " << setw(12) << dist << "  " << 0.0 << "\n";
//
//  Case 5, nonparallel, could intersect, and does not intersect.
//
//  L1: (1,1,1) + t * (0,1,2)
//  L2: (0,2,3) + t * (1,0,0)
//  lines intersect at (1,2,3), line segments do not.
//  Distance is 1.0;
//
  p1[0] = 1.0;
  p1[1] = 1.0;
  p1[2] = 1.0;

  p2[0] = 1.0;
  p2[1] = 4.0;
  p2[2] = 7.0;

  p3[0] =  0.0;
  p3[1] =  2.0;
  p3[2] =  3.0;

  p4[0] = -5.0;
  p4[1] =  2.0;
  p4[2] =  3.0;

  dist = segments_dist_3d ( p1, p2, p3, p4 );

  cout << "  " << setw(3) << 5
       << "  " << setw(12) << dist << "  " << 1.0 << "\n";
//
//  Case 6, nonparallel, can not intersect, "end-to-end".
//
//  L1: (2,2,1) + t * (0,1,2)  0 <= t <= 5
//  L2: (0,0,0) + t * (-1,-1,-1) 0 <= t <= 5
//  Distance is 3.
//
  p1[0] = 2.0;
  p1[1] = 2.0;
  p1[2] = 1.0;

  p2[0] =  2.0;
  p2[1] =  7.0;
  p2[2] = 11.0;

  p3[0] =  0.0;
  p3[1] =  0.0;
  p3[2] =  0.0;

  p4[0] =  -5.0;
  p4[1] =  -5.0;
  p4[2] =  -5.0;

  dist = segments_dist_3d ( p1, p2, p3, p4 );

  cout << "  " << setw(3) << 6
       << "  " << setw(12) << dist << "  " << 3.0 << "\n";
//
//  Case 7, nonparallel, can not intersect, "end-to-mid".
//
//  L1: (1,1,1) + t * (0,1,2) 0 <= t <= 5
//  L2: (0,4,7) + t * (-1,0,0) 0 <= t <= 5
//  Distance is 1.
//
  p1[0] = 1.0;
  p1[1] = 1.0;
  p1[2] = 1.0;

  p2[0] =  1.0;
  p2[1] =  6.0;
  p2[2] = 11.0;

  p3[0] =  0.0;
  p3[1] =  4.0;
  p3[2] =  7.0;

  p4[0] = -5.0;
  p4[1] =  4.0;
  p4[2] =  7.0;

  dist = segments_dist_3d ( p1, p2, p3, p4 );

  cout << "  " << setw(3) << 7
       << "  " << setw(12) << dist << "  " << 1.0 << "\n";
//
//  Case 8, nonparallel, can not intersect, "mid-to-mid".
//
//  L1: (0,5,10) + t * (1,-1,0) 0 <= t <= 5
//  L2: (0,0,0) + t * (1,1,0) 0 <= t <= 6
//  Distance = 10.
//
  p1[0] = 0.0;
  p1[1] = 5.0;
  p1[2] = 10.0;

  p2[0] =  5.0;
  p2[1] =  0.0;
  p2[2] =  10.0;

  p3[0] = 0.0;
  p3[1] = 0.0;
  p3[2] = 0.0;

  p4[0] = 6.0;
  p4[1] = 6.0;
  p4[2] = 0.0;

  dist = segments_dist_3d ( p1, p2, p3, p4 );

  cout << "  " << setw(3) << 8
       << "  " << setw(12) << dist << "  " << 10.0 << "\n";
//
//  Case 9, nonparallel, can not intersect, "mid-to-end".
//
//  L1: (-2,0,0) + t * (1,0,0) 0 <= t <= 12
//  L2: (-2,8,1) + t * (9,-4,-1) 0 <= t <= 1
//  Distance = 4.
//
  p1[0] = -2.0;
  p1[1] = 0.0;
  p1[2] = 0.0;

  p2[0] = 10.0;
  p2[1] =  0.0;
  p2[2] =  0.0;

  p3[0] = -2.0;
  p3[1] = 8.0;
  p3[2] = 1.0;

  p4[0] = 7.0;
  p4[1] = 4.0;
  p4[2] = 0.0;

  dist = segments_dist_3d ( p1, p2, p3, p4 );

  cout << "  " << setw(3) << 9
       << "  " << setw(12) << dist << "  " << 4.0 << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test044 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST044 tests SEGMENTS_INT_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 7

  double dist;
  int test;
  double p1;
  double p2;
  double q1;
  double q1_test[TEST_NUM] = { -1.0, 3.0, 1.0, 0.5, 0.25, 0.5, 2.0 };
  double q2;
  double q2_test[TEST_NUM] = { 1.0, 2.0, 2.0, -3.0, 0.50, 0.5, 2.0 };
  double r1;
  double r2;

  cout << "\n";
  cout << "TEST044\n";
  cout << "  SEGMENTS_INT_1D determines the intersection [R1,R2]\n";
  cout << "  of line segments [P1,P2] and [Q1,Q2] in 1D.\n";
  cout << "\n";
  cout << "  DIST is negative for overlap,\n";
  cout << "  0 for point intersection,\n";
  cout << "  positive if there is no overlap.\n";
  cout << "\n";
  cout << "  Test        P1        P2        Q1        Q2        R1        R2       DIST\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = -1.0;
    p2 = 1.0;
    q1 = q1_test[test];
    q2 = q2_test[test];

    dist = segments_int_1d ( p1, p2, q1, q2, &r1, &r2 );

    cout << "  " << setw(2) << test
         << "  " << setw(8) << p1
         << "  " << setw(8) << p2
         << "  " << setw(8) << q1
         << "  " << setw(8) << q2
         << "  " << setw(8) << r1
         << "  " << setw(8) << r2
         << "  " << setw(8) << dist << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test045 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST045 tests SEGMENTS_INT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  int flag;
  int i;
  double p1[DIM_NUM] = { -1.0, 3.0 };
  double p2[DIM_NUM] = {  1.0, 1.0 };
  double *q1;
  double q1_test[DIM_NUM*TEST_NUM] = {
    -1.0,  1.0,
     3.0, -1.0,
     0.0,  0.0,
     1.0,  2.0 };
  double *q2;
  double q2_test[DIM_NUM*TEST_NUM] = {
    1.0, -1.0,
    2.0,  0.0,
    0.0,  9.0,
    3.0,  2.0 };
  double r[DIM_NUM];
  int test;

  cout << "\n";
  cout << "TEST045\n";
  cout << "  SEGMENTS_INT_2D searches for an intersection of two\n";
  cout << "  line segments in 2D.\n";
  cout << "\n";
  cout << "  All tests use the same line segment 1:\n";
  cout << "\n";
  cout << "  P1:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p1[i];
  }
  cout << "\n";
  cout << "  P2:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p2[i];
  }

  for ( test = 0; test < TEST_NUM; test++ )
  {
    q1 = q1_test + test * DIM_NUM;
    q2 = q2_test + test * DIM_NUM;

    cout << "\n";
    cout << "\n";
    cout << "  Q1:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q1[i];
    }
    cout << "\n";
    cout << "  Q2:";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << q2[i];
    }
    cout << "\n";

    segments_int_2d ( p1, p2, q1, q2, &flag, r );

    if ( flag == 0 )
    {
      cout << "\n";
      cout << "  The line segments do not intersect.\n";
    }
    else if ( flag == 1 )
    {
      cout << "\n";
      cout << "  The line segments intersect at:\n";
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << "  " << setw(12) << r[i];
      }
      cout << "\n";
    }
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test046 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST046 tests MINABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2005
//
//  Author:
//
//    John Burkardt
//
{
  double x1;
  double x2;
  double x3;
  double xmin;
  double y1;
  double y2;
  double y3;
  double ymin;

  cout << "\n";
  cout << "TEST046\n";
  cout << "  MINABS finds the minimum of a function\n";
  cout << "    F(X) = a * ABS ( X ) + B\n";
  cout << "  within an interval, given three data points.\n";
//
//  Case 1: the three points lie on a straight line.
//  (XMIN=9,YMIN=2).
//
  x1 = 14.0;
  y1 = 7.0;

  x2 = 9.0;
  y2 = 2.0;

  x3 = 12.0;
  y3 = 5.0;

  minabs ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  cout << "\n";
  cout << "  The points lie on a straight line.\n";
  cout << "  XMIN = " << xmin << "  YMIN = " << ymin << "\n";
//
//  Case 2: the three points straddle a minimum.
//  (XMIN=7, YMIN=2).
//
  x1 = 3.0;
  y1 = 6.0;

  x2 = 12.0;
  y2 = 7.0;

  x3 = 9.0;
  y3 = 4.0;

  minabs ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  cout << "\n";
  cout << "  The points straddle a minimum.\n";
  cout << "  XMIN = " << xmin << "  YMIN = " << ymin << "\n";
//
//  Case 3: the three points straddle a maximum.
//  (XMIN=2, YMIN=5).
//
  x1 = 11.0;
  y1 = 6.0;

  x2 = 6.0;
  y2 = 9.0;

  x3 = 2.0;
  y3 = 5.0;

  minabs ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  cout << "\n";
  cout << "  The points straddle a maximum.\n";
  cout << "  XMIN = " << xmin << "  YMIN = " << ymin << "\n";

  return;
}
//****************************************************************************80

void test047 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST047 tests MINQUAD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  bool result;
  double x1;
  double x2;
  double x3;
  double xmin;
  double y1;
  double y2;
  double y3;
  double ymin;

  cout << "\n";
  cout << "TEST047\n";
  cout << "  MINQUAD finds the minimum of a function\n";
  cout << "    F(X) = A * X**2 + B * X + C\n";
  cout << "  within an interval, given three data points.\n";
//
//  Case 1: a minimum is in the interval.
//  y = ( x - 1 )**2 + 4
//
  x1 = 0.0;
  y1 = ( x1 - 1.0 ) * ( x1 - 1.0 ) + 4.0;

  x2 = 2.0;
  y2 = ( x2 - 1.0 ) * ( x2 - 1.0 ) + 4.0;

  x3 = 3.0;
  y3 = ( x3 - 1.0 ) * ( x3 - 1.0 ) + 4.0;

  result = minquad ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  if ( !result )
  {
    cout << "\n";
    cout << "Warning\n";
    cout << "  MINQUAD returned an error code.\n";
  }
  else
  {
    cout << "\n";
    cout << "  The minimum lies in the interval.\n";
    cout << "  X1, Y1 = " << x1 << "  " << y1 << "\n";
    cout << "  X2, Y2 = " << x2 << "  " << y2 << "\n";
    cout << "  X3, Y3 = " << x3 << "  " << y3 << "\n";
    cout << "  XMIN = " << xmin << " YMIN = " << ymin << "\n";
  }
//
//  Case 2: the minimum is to the left of the interval.
//  y = ( x - 1 )**2 + 4
//
  x1 = 2.0;
  y1 = ( x1 - 1.0 ) * ( x1 - 1.0 ) + 4.0;

  x2 = 4.0;
  y2 = ( x2 - 1.0 ) * ( x2 - 1.0 ) + 4.0;

  x3 = 5.0;
  y3 = ( x3 - 1.0 ) * ( x3 - 1.0 ) + 4.0;

  result = minquad ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  if ( !result )
  {
    cout << "\n";
    cout << "Warning\n";
    cout << "  MINQUAD returned an error code.\n";
  }
  else
  {
    cout << "\n";
    cout << "  The minimum is to the left of the interval\n";
    cout << "  X1, Y1 = " << x1 << "  " << y1 << "\n";
    cout << "  X2, Y2 = " << x2 << "  " << y2 << "\n";
    cout << "  X3, Y3 = " << x3 << "  " << y3 << "\n";
    cout << "  XMIN = " << xmin << " YMIN = " << ymin << "\n";
  }
//
//  Case 3: the function is flat.
//
  x1 = 11.0;
  y1 = 6.0;

  x2 = 6.0;
  y2 = 6.0;

  x3 = 2.0;
  y3 = 6.0;

  result = minquad ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  if ( !result )
  {
    cout << "\n";
    cout << "Warning\n";
    cout << "  MINQUAD returned an error code.\n";
  }
  else
  {
    cout << "\n";
    cout << "  The function is flat.\n";
    cout << "  X1, Y1 = " << x1 << "  " << y1 << "\n";
    cout << "  X2, Y2 = " << x2 << "  " << y2 << "\n";
    cout << "  X3, Y3 = " << x3 << "  " << y3 << "\n";
    cout << "  XMIN = " << xmin << " YMIN = " << ymin << "\n";
  }
//
//  Case 4: the function has a maximum.
//  y = - ( x - 1 )**2 + 4
//
  x1 = 0.0;
  y1 = - ( x1 - 1.0 ) * ( x1 - 1.0 ) + 4.0;

  x2 = 2.0;
  y2 = - ( x2 - 1.0 ) * ( x2 - 1.0 ) + 4.0;

  x3 = 3.0;
  y3 = - ( x3 - 1.0 ) * ( x3 - 1.0 ) + 4.0;

  result = minquad ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  if ( !result )
  {
    cout << "\n";
    cout << "Warning\n";
    cout << "  MINQUAD returned an error code.\n";
  }
  else
  {
    cout << "\n";
    cout << "  The function has a maximum.\n";
    cout << "  X1, Y1 = " << x1 << "  " << y1 << "\n";
    cout << "  X2, Y2 = " << x2 << "  " << y2 << "\n";
    cout << "  X3, Y3 = " << x3 << "  " << y3 << "\n";
    cout << "  XMIN = " << xmin << " YMIN = " << ymin << "\n";
  }

  return;
}
//****************************************************************************80

void test0475 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0475 tests OCTAHEDRON_SIZE_3D, OCTAHEDRON_SHAPE_3D, SHAPE_PRINT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int edge_num;
  int face_num;
  int *face_order;
  int face_order_max;
  int *face_point;
  int point_num;
  double *point_coord;

  cout << "\n";
  cout << "TEST0475\n";
  cout << "  For the octahedron,\n";
  cout << "  OCTAHEDRON_SIZE_3D returns dimension information;\n";
  cout << "  OCTAHEDRON_SHAPE_3D returns face and order information.\n";
  cout << "  SHAPE_PRINT_3D prints this information.\n";
//
//  Get the data sizes.
//
  octahedron_size_3d ( &point_num, &edge_num, &face_num, &face_order_max );

  cout << "\n";
  cout << "    Number of vertices: " << point_num << "\n";
  cout << "    Number of edges   : " << edge_num << "\n";
  cout << "    Number of faces   : " << face_num << "\n";
  cout << "    Maximum face order: " << face_order_max << "\n";
//
//  Make room for the data.
//
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];
  point_coord = new double[DIM_NUM*point_num];
//
//  Get the data.
//
  octahedron_shape_3d ( point_num, face_num, face_order_max, point_coord,
    face_order, face_point );
//
//  Print the data.
//
  shape_print_3d ( point_num, face_num, face_order_max,
    point_coord, face_order, face_point );

  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0477 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0477 tests PARALLELOGRAM_AREA_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double p[2*4] = {
    2.0, 7.0, 
    5.0, 7.0, 
    6.0, 9.0, 
    3.0, 9.0 };

  cout << "\n";
  cout << "TEST0477\n";
  cout << "  PARALLELOGRAM_AREA_2D finds the area of a\n";
  cout << "  parallelogram in 2D.\n";

  r8mat_transpose_print ( 2, 4, p, "  Vertices:" );

  area = parallelogram_area_2d ( p );

  cout << "\n";
  cout << "  AREA = " << area << "\n";

  return;
}
//****************************************************************************80

void test0478 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0478 tests PARALLELOGRAM_AREA_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double p[3*4] = {
    1.0,       2.0,       3.0, 
    2.4142137, 3.4142137, 3.0, 
    1.7071068, 2.7071068, 4.0, 
    0.2928931, 0.2928931, 4.0 };

  cout << "\n";
  cout << "TEST0478\n";
  cout << "  PARALLELOGRAM_AREA_3D finds the area of a\n";
  cout << "  parallelogram in 3D.\n";

  r8mat_transpose_print ( 3, 4, p, "  Vertices:" );

  area = parallelogram_area_3d ( p );

  cout << "\n";
  cout << "  AREA = " << area << "\n";

  return;
}
//****************************************************************************80

void test048 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST048 tests PARALLELOGRAM_CONTAINS_POINT_2D.
//
//  Discussion:
//
//    The four points are In, Out, Out, and Out.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2005
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  int i;
  bool inside;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
    1.0,  0.5,
    2.0,  0.0,
    0.5, -0.1,
    0.1,  0.5 };
  double p1[DIM_NUM] = { 0.0, 0.0 };
  double p2[DIM_NUM] = { 1.0, 0.0 };
  double p3[DIM_NUM] = { 1.0, 1.0 };
  int test;

  cout << "\n";
  cout << "TEST048\n";
  cout << "  PARALLELOGRAM_CONTAINS_POINT_2D determines if a point \n";
  cout << "    is within a parallelogram in 2D.\n";
  cout << "\n";
  cout << "       P     Inside?\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    inside = parallelogram_contains_point_2d ( p1, p2, p3, p );

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p[i];
    }
    cout << "  " << inside << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0485 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0485 tests PARALLELOGRAM_CONTAINS_POINT_2D
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 51

  int i;
  int j;
  double p[DIM_NUM];
  double p1[DIM_NUM] = { 0.2, 0.0 };
  double p2[DIM_NUM] = { 0.4, 0.6 };
  double p3[DIM_NUM] = { 0.6, 0.4 };
  double xhi =  1.0;
  double xlo =  0.0;
  double yhi =  1.0;
  double ylo =  0.0;

  cout << "\n";
  cout << "TEST0485\n";
  cout << "  PARALLELOGRAM_CONTAINS_POINT_2D reports if a parallelogram\n";
  cout << "  contains a point.\n";
  cout << "\n";
  cout << "  We will call the function repeatedly, and draw\n";
  cout << "  a sketch of the unit square.\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    p[1] = ( ( N - i     ) * yhi
           + (     i - 1 ) * ylo )
           / ( N     - 1 );

    cout << "  ";
    for ( j = 1; j <= N; j++ )
    {
      p[0] = ( ( N - j     ) * xlo
             + (     j - 1 ) * xhi )
             / ( N     - 1 );

      if ( parallelogram_contains_point_2d ( p1, p2, p3, p ) )
      {
        cout << "*";
      }
      else
      {
        cout << "-";
      }
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test049 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST049 tests PARALLELOGRAM_CONTAINS_POINT_3D.
//
//  Discussion:
//
//    The points are In, Out, Out, Out, Out
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 5

  int i;
  bool inside;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
    1.0,  1.0,  0.5,
    3.0,  3.0,  0.0,
    0.5,  0.5, -0.1,
    0.1,  0.1,  0.5,
    1.5,  1.6,  0.5 };
  double p1[DIM_NUM] = { 0.0, 0.0, 0.0 };
  double p2[DIM_NUM] = { 2.0, 2.0, 0.0 };
  double p3[DIM_NUM] = { 1.0, 1.0, 1.0 };
  int test;

  cout << "\n";
  cout << "TEST049\n";
  cout << "  PARALLELOGRAM_CONTAINS_POINT_3D determines if a point\n";
  cout << "    is within a parallelogram in 3D.\n";
  cout << "\n";
  cout << "           P           Inside?\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    inside = parallelogram_contains_point_3d ( p1, p2, p3, p );

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p[i];
    }
    cout << "  " << inside << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0493 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0493 tests PARABOLA_EX and PARABOLA_EX2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  int ierror = 0;
  double x1;
  double x2;
  double x3;
  double xmin;
  double y1;
  double y2;
  double y3;
  double ymin;

  cout << "\n";
  cout << "TEST0493\n";
  cout << "  PARABOLA_EX finds the extreme value of a parabola\n";
  cout << "    determined by three points.\n";
  cout << "  PARABOLA_EX2 finds the extreme value of a parabola\n";
  cout << "    determined by three points.\n";

  a =  2.0;
  b = -4.0;
  c = 10.0;

  x1 = 1.0;
  y1 = a * x1 * x1 + b * x1 + c;
  x2 = 2.0;
  y2 = a * x2 * x2 + b * x2 + c;
  x3 = 3.0;
  y3 = a * x3 * x3 + b * x3 + c;

  cout << "\n";
  cout << "  Parabolic coefficients (A,B,C) = "
    << a << "  " << b << "  " << c << "\n";
  cout << "\n";
  cout << "  X, Y data:\n";
  cout << "\n";
  cout << "  X1, Y1 = " << x1 << "  " << y1 << "\n";
  cout << "  X2, Y2 = " << x2 << "  " << y2 << "\n";
  cout << "  X3, Y3 = " << x3 << "  " << y3 << "\n";

  a = 0.0;
  b = 0.0;
  c = 0.0;

  ierror = parabola_ex ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  if ( ierror == 0 )
  {
    cout << "\n";
    cout << "  PARABOLA_EX returns (XMIN,YMIN) = "
         << xmin << "  " << ymin << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  PARABOLA_EX returns error code " << ierror << "\n";
  }

  ierror = parabola_ex2 ( x1, y1, x2, y2, x3, y3, &xmin, &ymin, &a, &b, &c );

  if ( ierror == 0 )
  {
    cout << "\n";
    cout << "  PARABOLA_EX2 returns (XMIN,YMIN) = " << xmin
      << "  " << ymin << "\n";
    cout << "  and (A,B,C) = " << a << "  " << b << "  " << c << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  PARABOLA_EX2 returns error code " << ierror << "\n";
  }

  return;
}
//****************************************************************************80

void test0495 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0495 tests PARALLELEPIPED_POINT_DIST_3D.
//
//  Discussion:
//
//    The points tested are:
//
//    1: Center of box.
//    2: The middle of a face.
//    3: The middle of an edge.
//    4: A corner.
//    5: Close to a face.
//    6: Close to an edge.
//    7: Close to a corner.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 7

  double dist;
  int i;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
    1.0,  4.0,  0.5,
    1.0,  0.0,  0.5,
    0.0,  4.0,  1.0,
    2.0,  8.0,  1.0,
   -0.5,  4.0,  0.5,
    1.0, -1.0, -1.0,
    3.0,  9.0,  2.0 };
  double p1[DIM_NUM] = { 0.0, 0.0, 0.0 };
  double p2[DIM_NUM] = { 2.0, 0.0, 0.0 };
  double p3[DIM_NUM] = { 0.0, 8.0, 0.0 };
  double p4[DIM_NUM] = { 0.0, 0.0, 1.0 };
  int test;

  cout << "\n";
  cout << "TEST0495\n";
  cout << "  PARALLELEPIPED_POINT_DIST_3D computes the distance\n";
  cout << "    from a point to a box (parallelipiped) in 3D.\n";
  cout << "\n";
  cout << "  The 4 box corners that are specified:\n";
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p1[i];
  }
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p2[i];
  }
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p3[i];
  }
  cout << "\n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p4[i];
  }
  cout << "\n";
  cout << "\n";
  cout << " TEST          P                  Distance to box\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    dist = parallelepiped_point_dist_3d ( p1, p2, p3, p4, p );

    cout << "  " << setw(3) << test;
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << p[i];
    }
    cout << "  " << setw(12) << dist << "\n";

  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test050 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST050 tests PLANE_EXP_NORMAL_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double p1[DIM_NUM] = { -10.56, -10.56, 78.09 };
  double p2[DIM_NUM] = {  44.66, -65.77,   0.0 };
  double p3[DIM_NUM] = {  44.66,  44.66,   0.0 };
  double pn[DIM_NUM];

  cout << "\n";
  cout << "TEST050\n";
  cout << "  PLANE_EXP_NORMAL_3D finds the normal to a plane.\n";
  cout << "\n";
  cout << "  Coordinates of 3 points:\n";
  cout << "\n";
  cout << "  P1 = " << setw(12) << p1[0]
       << "  "      << setw(12) << p1[1]
       << "  "      << setw(12) << p1[2] << "\n";
  cout << "  P2 = " << setw(12) << p2[0]
       << "  "      << setw(12) << p2[1]
       << "  "      << setw(12) << p2[2] << "\n";
  cout << "  P3 = " << setw(12) << p3[0]
       << "  "      << setw(12) << p3[1]
       << "  "      << setw(12) << p3[2] << "\n";

  plane_exp_normal_3d ( p1, p2, p3, pn );

  cout << "\n";
  cout << "  Unit normal vector:\n";
  cout << "\n";
  cout << "  PN = " << setw(12) << pn[0]
       << "  "      << setw(12) << pn[1]
       << "  "      << setw(12) << pn[2] << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test051 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST051 tests PLANE_EXP2IMP_3D.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double a;
  double b;
  double c;
  double d;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];

  cout << "\n";
  cout << "TEST051\n";
  cout << "  PLANE_EXP2IMP_3D puts a plane defined by\n";
  cout << "    3 points into A*X+B*Y+C*Z+D = 0 form.\n";
  cout << "\n";

  p1[0] = -1.0;
  p1[1] = 0.0;
  p1[2] = -1.0;

  p2[0] = -4.0;
  p2[1] = 0.0;
  p2[2] = 0.0;

  p3[0] = -20.0;
  p3[1] = 2.0;
  p3[2] = 4.0;

  plane_exp2imp_3d ( p1, p2, p3, &a, &b, &c, &d );

  cout << "\n";
  cout << "  Input:\n";
  cout << "\n";
  cout << "  P1 = " << setw(12) << p1[0]
       << "  "      << setw(12) << p1[1]
       << "  "      << setw(12) << p1[2] << "\n";
  cout << "  P2 = " << setw(12) << p2[0]
       << "  "      << setw(12) << p2[1]
       << "  "      << setw(12) << p2[2] << "\n";
  cout << "  P3 = " << setw(12) << p3[0]
       << "  "      << setw(12) << p3[1]
       << "  "      << setw(12) << p3[2] << "\n";
  cout << "\n";
  cout << "  Output:\n";
  cout << "\n";
  cout << "  (A,B,C,D)= " << a << "  " << b << "  " << c << "  " << d << "\n";
  cout << "  Correct answer is a multiple of 1, 2, 3, 4.\n";

  p1[0] = -16.0;
  p1[1] = 2.0;
  p1[2] = 4.0;

  p2[0] = 0.0;
  p2[1] = 0.0;
  p2[2] = 0.0;

  p3[0] = 4.0;
  p3[1] = -2.0;
  p3[2] = 0.0;

  plane_exp2imp_3d ( p1, p2, p3, &a, &b, &c, &d );

  cout << "\n";
  cout << "  Input:\n";
  cout << "\n";
  cout << "  P1 = " << setw(12) << p1[0]
       << "  "      << setw(12) << p1[1]
       << "  "      << setw(12) << p1[2] << "\n";
  cout << "  P2 = " << setw(12) << p2[0]
       << "  "      << setw(12) << p2[1]
       << "  "      << setw(12) << p2[2] << "\n";
  cout << "  P3 = " << setw(12) << p3[0]
       << "  "      << setw(12) << p3[1]
       << "  "      << setw(12) << p3[2] << "\n";
  cout << "\n";
  cout << "  Output:\n";
  cout << "\n";
  cout << "  (A,B,C,D)= " << a << "  " << b << "  " << c << "  " << d << "\n";
  cout << "  Correct answer is a multiple of 1, 2, 3, 0.\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test052 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST052 tests PLANE_EXP2NORMAL_3D.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int i;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];
  double pn[DIM_NUM];
  double pp[DIM_NUM];

  cout << "\n";
  cout << "TEST052\n";
  cout << "  PLANE_EXP2NORMAL_3D puts a plane defined by\n";
  cout << "    3 points into point, normal form\n";
  cout << "\n";

  p1[0] = -1.0;
  p1[1] = 0.0;
  p1[2] = -1.0;

  p2[0] = -4.0;
  p2[1] = 0.0;
  p2[2] = 0.0;

  p3[0] = -20.0;
  p3[1] = 2.0;
  p3[2] = 4.0;

  plane_exp2normal_3d ( p1, p2, p3, pp, pn );

  cout << "\n";
  cout << "  Input:\n";
  cout << "\n";
  cout << "  P1 = " << setw(12) << p1[0]
       << "  "      << setw(12) << p1[1]
       << "  "      << setw(12) << p1[2] << "\n";
  cout << "  P2 = " << setw(12) << p2[0]
       << "  "      << setw(12) << p2[1]
       << "  "      << setw(12) << p2[2] << "\n";
  cout << "  P3 = " << setw(12) << p3[0]
       << "  "      << setw(12) << p3[1]
       << "  "      << setw(12) << p3[2] << "\n";
  cout << "\n";
  cout << "  Output:\n";
  cout << "\n";
  cout << "  PP = " << pp[0] << "  " << pp[1] << "  " << pp[2] << "\n";
  cout << "  PN = " << pn[0] << "  " << pn[1] << "  " << pn[2] << "\n";

  p1[0] = -16.0;
  p1[1] = 2.0;
  p1[2] = 4.0;

  p2[0] = 0.0;
  p2[1] = 0.0;
  p2[2] = 0.0;

  p3[0] = 4.0;
  p3[1] = -2.0;
  p3[2] = 0.0;

  cout << "  P1: ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p1[i];
  }
  cout << "\n";
  cout << "  P2: ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p2[i];
  }
  cout << "\n";
  cout << "  P3: ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p3[i];
  }
  cout << "\n";

  plane_exp2normal_3d ( p1, p2, p3, pp, pn );

  cout << "\n";
  cout << "  Output:\n";
  cout << "\n";
  cout << "  PP = " << pp[0] << "  " << pp[1] << "  " << pp[2] << "\n";
  cout << "  PN = " << pn[0] << "  " << pn[1] << "  " << pn[2] << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test053 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST053 tests PLANE_EXP_PROJECT_3D.
//
//  Discussion:
//
//    1: Projection   is ( 0, 0.5, 0.5 ), IVIS is 3.
//    2: Projection   is ( 4, 5, -8 ), IVIS is 2.
//    3: Projection   is ( 0.33, 0.33, 0.33), IVIS is 1.
//    4: "Projection" is ( 0, 0, 0 ), IVIS is 0.
//    5: Projection   is ( 1, 0, 0 ), IVIS is -1.
//
//  Modified:
//
//    18 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define N 5

  int i;
  int ivis[N];
  int j;
  double p1[DIM_NUM] = { 1.0, 0.0, 0.0 };
  double p2[DIM_NUM] = { 0.0, 1.0, 0.0 };
  double p3[DIM_NUM] = { 0.0, 0.0, 1.0 };
  double pf[DIM_NUM] = { 0.0, 0.0, 0.0 };
  double po[DIM_NUM*N] = {
     0.00,  2.00,  2.00,
     4.00,  5.00, -8.00,
     0.25,  0.25,  0.25,
     5.00, -2.00, -3.00,
    -2.00,  0.00,  0.00 };
  double pp[DIM_NUM*N];

  cout << "\n";
  cout << "TEST053\n";
  cout << "  PLANE_EXP_PROJECT_3D projects a point through\n";
  cout << "    a focus point into a plane.\n";
  cout << "\n";
  cout << "      PO      PP      IVIS\n";
  cout << "\n";

  plane_exp_project_3d ( p1, p2, p3, pf, N, po, pp, ivis );

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(10) << po[i+j*DIM_NUM];
    }
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(10) << pp[i+j*DIM_NUM];
    }
    cout << "  " << setw(4) << ivis[j];
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test054 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST054 tests PLANE_IMP2EXP_3D.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double a;
  double b;
  double c;
  double d;
  int i;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];

  cout << "\n";
  cout << "TEST054\n";
  cout << "  PLANE_IMP2EXP_3D converts a plane in implicit\n";
  cout << "    (A,B,C,D) form to explicit form.\n";

  a = 1.0;
  b = -2.0;
  c = -3.0;
  d = 6.0;

  cout << "\n";
  cout << "  A = " << a
       << "  B = " << b
       << "  C = " << c
       << "  D = " << d << "\n";

  plane_imp2exp_3d ( a, b, c, d, p1, p2, p3 );

  cout << "  P1: ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p1[i];
  }
  cout << "\n";
  cout << "  P2: ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p2[i];
  }
  cout << "\n";
  cout << "  P3: ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p3[i];
  }
  cout << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test055 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST055 tests PLANE_IMP2NORMAL_3D.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double a;
  double b;
  double c;
  double d;
  double pn[DIM_NUM];
  double pp[DIM_NUM];

  cout << "\n";
  cout << "TEST055\n";
  cout << "  PLANE_IMP2NORMAL_3D converts a plane in implicit\n";
  cout << "    (A,B,C,D) form to point, normal form.\n";

  a = 1.0;
  b = -2.0;
  c = -3.0;
  d = 6.0;

  cout << "\n";
  cout << "  A = " << a
       << "  B = " << b
       << "  C = " << c
       << "  D = " << d << "\n";

  plane_imp2normal_3d ( a, b, c, d, pp, pn );

  cout << "\n";
  cout << "  PP = " << pp[0] << "  " << pp[1] << "  " << pp[2] << "\n";
  cout << "  PN = " << pn[0] << "  " << pn[1] << "  " << pn[2] << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test056 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST056 tests PLANE_IMP_LINE_PAR_INT_3D.
//
//  Modified:
//
//    18 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double a;
  double b;
  double c;
  double d;
  double f;
  double g;
  double h;
  int i;
  bool intersect;
  double p[DIM_NUM];
  double x0;
  double y0;
  double z0;

  cout << "\n";
  cout << "TEST056\n";
  cout << "  PLANE_IMP_LINE_PAR_INT_3D finds the\n";
  cout << "    intersection of an implicit plane and\n";
  cout << "    a parametric line, in 3D.\n";

  a = 1.0;
  b = -2.0;
  c = -3.0;
  d = 6.0;

  f = 2.0;
  g = 1.0;
  h = 5.0;
  x0 = 3.0;
  y0 = 0.0;
  z0 = -7.0;

  intersect = plane_imp_line_par_int_3d ( a, b, c, d, x0, y0, z0, f, g, h, p );

  if ( intersect )
  {
    cout << "\n";
    cout << "  The plane and line intersect at\n";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p[i];
    }
    cout << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  The plane and the line do not intersect.\n";
  }

  cout << "\n";
  cout << "  Expected answer:\n";
  cout << "    The plane and line intersect at \n";
  cout << "    7, 2, 3.\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test057 ( )

//****************************************************************************80*
//
//  Purpose:
//
//    TEST057 tests PLANE_IMP_SEGMENT_NEAR_3D.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  double a;
  double b;
  double c;
  double d;
  double dist;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double pls[DIM_NUM];
  double pp[DIM_NUM];
  int test;

  cout << "\n";
  cout << "TEST057\n";
  cout << "  PLANE_IMP_SEGMENT_NEAR_3D finds the point\n";
  cout << "    on a line segment nearest a plane.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1[0] = 3.0;
    p1[1] = 0.0;
    p1[2] = -7.0;

    if ( test == 0 )
    {
      p2[0] = 9.0;
      p2[1] = 3.0;
      p2[2] = 8.0;
    }
    else if ( test == 1 )
    {
      p2[0] = 5.0;
      p2[1]= 1.0;
      p2[2] = -2.0;
    }

    a = 1.0;
    b = -2.0;
    c = -3.0;
    d = 6.0;

    plane_imp_segment_near_3d ( p1, p2, a, b, c, d, &dist, pp, pls );

    cout << "\n";
    cout << "  The distance between the plane and the\n";
    cout << "  line segment is " << dist << "\n";
    cout << "\n";
    cout << "  A nearest point on the line segment is\n";
    cout << "  PLS = " << pls[0] << "  " << pls[1] << "  " << pls[2] << "\n";
    cout << "\n";
    cout << "  A nearest point on the plane is\n";
    cout << "  PP = " << pp[0] << "  " << pp[1] << "  " << pp[2] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test058 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST058 tests PLANE_IMP_POINT_DIST_3D, PLANE_IMP_POINT_DIST_SIGNED_3D;
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 4

  double a;
  double b;
  double c;
  double d;
  double dist;
  double dist_signed;
  double p[3];
  int test;
  double xtest[TEST_NUM];
  double ytest[TEST_NUM];
  double ztest[TEST_NUM];

  a =    0.0;
  b =    0.0;
  c =    1.0;
  d = - 10.0;

  xtest[0] = - 12.0;
  ytest[0] =   14.0;
  ztest[0] =    0.0;

  xtest[1] =    7.0;
  ytest[1] =    8.0;
  ztest[1] =    9.0;

  xtest[2] =    1.0;
  ytest[2] =    2.0;
  ztest[2] =   10.0;

  xtest[3] =    0.0;
  ytest[3] =    0.0;
  ztest[3] =   12.0;

  cout << "\n";
  cout << "TEST058\n";
  cout << "  PLANE_IMP_POINT_DIST_3D computes the distance\n";
  cout << "    between an implicit plane and a point in 3D;\n";
  cout << "  PLANE_IMP_POINT_DIST_SIGNED 3D computes the\n";
  cout << "    signed distance between an implicit plane\n";
  cout << "    and a point in 3D.\n";
  cout << "\n";
  cout << "  For all tests, we use the implicit plane with\n";
  cout << "  (A,B,C,D) = "
    << a << "  " << b << "  " << c << "  " << d << "\n";
  cout << "\n";
  cout << "  P     DISTANCE   SIGNED_DISTANCE\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p[0] = xtest[test];
    p[1] = ytest[test];
    p[2] = ztest[test];

    dist = plane_imp_point_dist_3d ( a, b, c, d, p );

    dist_signed = plane_imp_point_dist_signed_3d ( a, b, c, d, p );

    cout << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << p[2]
         << "  " << setw(10) << dist
         << "  " << setw(10) << dist_signed << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test059 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST059 tests PLANE_IMP_TRIANGLE_NEAR_3D.
//
//  Modified:
//
//    18 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  double a;
  double b;
  double c;
  double d;
  double dist;
  int near_num;
  double pn[DIM_NUM*6];
  double t[DIM_NUM*3] = {
     3.0,  0.0, -7.0,
    13.0, -4.0, -1.0,
     0.0,  0.0,  0.0 };
  int test;

  cout << "\n";
  cout << "TEST059\n";
  cout << "  PLANE_IMP_TRIANGLE_NEAR_3D finds the nearest\n";
  cout << "    points on an implicit plane and a triangle.\n";

  a = 1.0;
  b = -2.0;
  c = -3.0;
  d = 6.0;

  cout << "\n";
  cout << "  Implicit plane: A*X + B*Y + C*Z + D = 0.\n";
  cout << "  A = " << a
       << "  B = " << b
       << "  C = " << c
       << "  D = " << d << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    if ( test == 1 )
    {
      t[0+2*DIM_NUM] =  5.0;
      t[1+2*DIM_NUM] =  1.0;
      t[2+2*DIM_NUM] = -2.0;
    }
    else if ( test == 2 )
    {
      t[0+2*DIM_NUM] = 9.0;
      t[1+2*DIM_NUM] = 3.0;
      t[2+2*DIM_NUM] = 8.0;
    }

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    near_num = plane_imp_triangle_near_3d ( t, a, b, c, d, &dist, pn );

    cout << "\n";
    cout << "  Triangle to plane distance is " << dist << "\n";

    r8mat_transpose_print ( DIM_NUM, near_num, pn, "  Nearest points:" );
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test060 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST060 tests PLANE_IMP_TRIANGLE_INT_3D.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 4

  double a;
  double b;
  double c;
  double d;
  int int_num;
  int j;
  double p[DIM_NUM*3];
  double t[DIM_NUM*3];
  int test;

  cout << "\n";
  cout << "TEST060\n";
  cout << "  PLANE_IMP_TRIANGLE_INT_3D finds the\n";
  cout << "    intersection points of an implicit plane\n";
  cout << "    and a triangle.\n";

  a = 1.0;
  b = -2.0;
  c = -3.0;
  d = 6.0;

  cout << "\n";
  cout << "  The implicit plane: A*X + B*Y + C*Z + D = 0.\n";
  cout << "  (A,B,C,D) = " << a << "  " << b << "  " << c << "  " << d << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    if ( test == 0 )
    {
      t[0+0*3] =  3.0;
      t[1+0*3] =  0.0;
      t[2+0*3] = -7.0;
      t[0+1*3] = 13.0;
      t[1+1*3] = -4.0;
      t[2+1*3] = -1.0;
      t[0+2*3] =  5.0;
      t[1+2*3] =  1.0;
      t[2+2*3] = -2.0;
    }
    else if ( test == 1 )
    {
      t[0+0*3] =  3.0;
      t[1+0*3] =  0.0;
      t[2+0*3] = -7.0;
      t[0+1*3] = 13.0;
      t[1+1*3] = -4.0;
      t[2+1*3] = -1.0;
      t[0+2*3] =  9.0;
      t[1+2*3] =  3.0;
      t[2+2*3] =  8.0;
    }
    else if ( test == 2 )
    {
      t[0+0*3] = -6.0;
      t[1+0*3] =  0.0;
      t[2+0*3] =  0.0;
      t[0+1*3] =  0.0;
      t[1+1*3] =  3.0;
      t[2+1*3] =  0.0;
      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] =  2.0;
    }
    else if ( test == 3 )
    {
      t[0+0*3] = -4.0;
      t[1+0*3] = +1.0;
      t[2+0*3] =  0.0;
      t[0+1*3] =  0.0;
      t[1+1*3] =  6.0;
      t[2+1*3] = -2.0;
      t[0+2*3] =  0.0;
      t[1+2*3] =  0.0;
      t[2+2*3] =  1.0;
    }

    cout << "\n";
    cout << "  Case " << test << "\n";

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    plane_imp_triangle_int_3d ( a, b, c, d, t, &int_num, p );

    cout << "\n";
    cout << "  Number of intersection points is " << int_num << "\n";
    cout << "\n";

    if ( 0 < int_num )
    {
      for ( j = 0; j < int_num; j++ )
      {
        cout << "  " << setw(4)  << j+1
             << "  " << setw(10) << p[0+j*3]
             << "  " << setw(10) << p[1+j*3]
             << "  " << setw(10) << p[2+j*3] << "\n";
      }
    }

  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test061 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST061 tests PLANE_NORMAL_BASIS_3D.
//
//  Modified:
//
//    27 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 5

  double b[DIM_NUM*DIM_NUM];
  double *pn;
  double *pp;
  double pq[DIM_NUM];
  double pr[DIM_NUM];
  int seed = 123456789;
  int test;

  cout << "\n";
  cout << "TEST061\n";
  cout << "  PLANE_NORMAL_BASIS_3D, given a plane in\n";
  cout << "    point, normal form (P,N), finds two unit\n";
  cout << "    vectors Q and R that lie in the plane\n";
  cout << "    and are mutually orthogonal.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    pp = r8vec_uniform_01_new ( DIM_NUM, seed );
    pn = r8vec_uniform_01_new ( DIM_NUM, seed );

    plane_normal_basis_3d ( pp, pn, pq, pr );

    if ( test == 0 )
    {
      cout << "\n";
      cout << "  Data for first case:\n";
      cout << "\n";
      cout << "  PP = " << pp[0] << "  " << pp[1] << "  " << pp[2] << "\n";
      cout << "  PN = " << pn[0] << "  " << pn[1] << "  " << pn[2] << "\n";
      cout << "  PQ = " << pq[0] << "  " << pq[1] << "  " << pq[2] << "\n";
      cout << "  PR = " << pr[0] << "  " << pr[1] << "  " << pr[2] << "\n";
      cout << "\n";
    }
    b[0+0*3] = r8vec_dot_product ( DIM_NUM, pn, pn );
    b[1+0*3] = r8vec_dot_product ( DIM_NUM, pn, pq );
    b[2+0*3] = r8vec_dot_product ( DIM_NUM, pn, pr );
    b[0+1*3] = r8vec_dot_product ( DIM_NUM, pq, pn );
    b[1+1*3] = r8vec_dot_product ( DIM_NUM, pq, pq );
    b[2+1*3] = r8vec_dot_product ( DIM_NUM, pq, pr );
    b[0+2*3] = r8vec_dot_product ( DIM_NUM, pr, pn );
    b[1+2*3] = r8vec_dot_product ( DIM_NUM, pr, pq );
    b[2+2*3] = r8vec_dot_product ( DIM_NUM, pr, pr );

    r8mat_print ( DIM_NUM, DIM_NUM, b, "  Dot product matrix:" );

    delete [] pp;
    delete [] pn;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0615 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0615 tests PLANE_NORMAL_LINE_EXP_INT_3D.
//
//  Modified:
//
//    18 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int i;
  int ival;
  double normal[DIM_NUM] = { 1.0, -2.0, -3.0 };
  double p1[DIM_NUM] = { 3.0, 0.0, -7.0 };
  double p2[DIM_NUM] = { 5.0, 1.0, -2.0 };
  double pint[DIM_NUM];
  double pp[DIM_NUM] = { -1.0, +1.0, +1.0 };
  double temp;

  cout << "\n";
  cout << "TEST0615\n";
  cout << "  PLANE_NORMAL_LINE_EXP_INT_3D finds the\n";
  cout << "    intersection of a normal plane and\n";
  cout << "    an explicit line, in 3D.\n";

  temp = 0.0;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    temp = temp + pow ( normal[i], 2 );
  }
  temp = sqrt ( temp );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    normal[i] = normal[i] / temp;
  }

  r8vec_print ( DIM_NUM, pp, "  Plane point PP:" );
  r8vec_print ( DIM_NUM, normal, "  Plane Normal:" );

  r8vec_print ( DIM_NUM, p1, "  Line point P1:" );
  r8vec_print ( DIM_NUM, p2, "  Line point P2:" );

  ival = plane_normal_line_exp_int_3d ( pp, normal, p1, p2, pint );

  cout << "\n";

  if ( ival == 0 )
  {
    cout << "  The plane and line do not intersect.\n";
  }
  else if ( ival == 1 )
  {
    r8vec_print ( DIM_NUM, pint, "  The unique intersection point:" );
  }
  else if ( ival == 2 )
  {
    r8vec_print ( DIM_NUM, pint, "  One of infinitely many intersections:" );
  }
  else
  {
    cout << "  The plane and the line do not intersect.\n";
  }

  cout << "\n";
  cout << "  Expected answer:\n";
  cout << "    The unique intersection point:\n";
  cout << "    7, 2, 3.\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0616 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0616 tests PLANE_NORMAL_QR_TO_XYZ and PLANE_NORMAL_XYZ_TO_QR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  double dif;
  int i;
  int j;
  int m = 3;
  int n = 5;
  double *normal;
  double *pp;
  double pq[3];
  double pr[3];
  double *qr1;
  double *qr2;
  int seed;
  double t;
  double *xyz;

  seed = 123456789;

  cout << "\n";
  cout << "TEST0616\n";
  cout << "  For a normal plane, with point PP and NORMAL vector,\n";
  cout << "  and in-plane basis vectors PQ and PR,\n";
  cout << "  PLANE_NORMAL_QR_TO_XYZ converts QR to XYZ coordinates;\n";
  cout << "  PLANE_NORMAL_XYZ_TO_QR converts XYZ to QR coordinates.\n";
//
//  Choose PP and NORMAL at random.
//
  pp = r8vec_uniform_01_new ( m, seed );

  normal = r8vec_uniform_01_new ( m, seed );
//
//  Compute in-plane basis vectors PQ and PR.
//
  plane_normal_basis_3d ( pp, normal, pq, pr );
//
//  Choose random Q, R coordinates.
//
  qr1 = r8mat_uniform_01_new ( m - 1, n, seed );
//
//  Convert to XYZ.
//
  xyz = plane_normal_qr_to_xyz ( pp, normal, pq, pr, n, qr1 );
//
//  Convert XYZ to QR.
//
  qr2 = plane_normal_xyz_to_qr ( pp, normal, pq, pr, n, xyz );

  dif = 0.0;
  for ( j = 0; j < n; j++ )
  {
    t = 0.0;
    for ( i = 0; i < m - 1; i++ )
    {
      t = t + pow ( qr1[0+j*2] - qr2[0+j*2], 2 );
    }
    t = sqrt ( t );
    dif = r8_max ( dif, t );
  }

  cout << "\n";
  cout << "  Maximum difference was " << dif << "\n";

  delete [] normal;
  delete [] pp;
  delete [] qr1;
  delete [] qr2;
  delete [] xyz;

  return;
}
//****************************************************************************80

void test0617 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0617 tests PLANE_NORMAL_TETRAHEDRON_INTERSECT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  int l;
  int int_num;
  double normal[3];
  double pint[3*4];
  double pp[3];
  double t[3*4] = {
    0.0, 0.0, 0.0, 
    1.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 
    0.0, 0.0, 1.0 };

  cout << "\n";
  cout << "TEST0617\n";
  cout << "  PLANE_NORMAL_TETRAHEDRON_INTERSECT determines\n";
  cout << "  the intersection of a plane and tetrahedron.\n";

  for ( k = 1; k <= 2; k++ )
  {
    if ( k == 1 )
    {
      normal[0] = 0.0;
      normal[1] = 0.0;
      normal[2] = 1.0;
    }
    else
    {
      normal[0] = 1.0 / sqrt ( 2.0 );
      normal[1] = 1.0 / sqrt ( 2.0 );
      normal[2] = 0.0;
    }

    cout << "\n";
    cout << "  Plane normal vector number " << k << "\n";
    cout << "\n";
    for ( i = 0; i < 3; i++ )
    {
      cout << "  " << normal[i];
    }
    cout << "\n";

    for ( l = 0; l <= 6; l++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        pp[i] = normal[i] * ( double ) ( l ) / 5.0;
      }
      plane_normal_tetrahedron_intersect ( pp, normal, t, &int_num, pint );

      cout << "\n";
      cout << "  Point on plane:\n";
      cout << "\n";
      for ( i = 0; i < 3; i++ )
      {
        cout << "  " << pp[i];
      }
      cout << "\n";
      cout << "\n";
      cout << "  Number of intersection points = " << int_num << "\n";
      cout << "\n";
      for ( j = 0; j < int_num; j++ )
      {
        cout << "  " << setw(4) << j;
        for ( i = 0; i < 3; i++ )
        {
          cout << "  " << pint[i+j*3];
        }
        cout << "\n";
      }
    }
  }
  return;
}
//****************************************************************************80

void test062 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST062 tests PLANE_NORMAL_TRIANGLE_INT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 4

  int int_num;
  double normal[DIM_NUM] = { 1.0, -2.0, -3.0 };
  double pint[DIM_NUM*3];
  double pp[DIM_NUM] = { 0.0, 0.0, 2.0 };
  double *t;
  double t_test[DIM_NUM*3*TEST_NUM] = {
         3.0,  0.0, -7.0,
        13.0, -4.0, -1.0,
         5.0,  1.0, -2.0,
         3.0,  0.0, -7.0,
        13.0, -4.0, -1.0,
         9.0,  3.0,  8.0,
        -6.0,  0.0,  0.0,
         0.0,  3.0,  0.0,
         0.0,  0.0,  2.0,
        -4.0,  1.0,  0.0,
         0.0,  6.0, -2.0,
         0.0,  0.0,  1.0 };
  int test;

  cout << "\n";
  cout << "TEST062\n";
  cout << "  PLANE_NORMAL_TRIANGLE_INT_3D finds the\n";
  cout << "    intersection points of a normal form plane\n";
  cout << "    and a triangle.\n";

  r8vec_print ( DIM_NUM, pp, "  The point PP:" );
  r8vec_print ( DIM_NUM, normal, "  The normal vector N:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    t = t_test + DIM_NUM * 3 * test;

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    int_num = plane_normal_triangle_int_3d ( pp, normal, t, pint );

    cout << "\n";
    cout << "  Number of intersection points is " << int_num << "\n";
    cout << "\n";

    r8mat_transpose_print ( DIM_NUM, int_num, pint, "  Intersection points:" );
  }
  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test063 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST063 tests PLANE_NORMAL2EXP_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int i;
  double normal[DIM_NUM] = { -0.2672612, -0.5345225, -0.8017837 };
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];
  double pp[DIM_NUM] = { -1.0,   0.0, -1.0 };

  cout << "\n";
  cout << "TEST063\n";
  cout << "  PLANE_NORMAL2EXP_3D puts a plane defined by\n";
  cout << "    point, normal form into explicit form.\n";

  r8vec_print ( DIM_NUM, pp, "  The point PP:" );

  r8vec_print ( DIM_NUM, normal, "  Normal vector:" );

  plane_normal2exp_3d ( pp, normal, p1, p2, p3 );

  cout << "  P1: ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p1[i];
  }
  cout << "\n";
  cout << "  P2: ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p2[i];
  }
  cout << "\n";
  cout << "  P3: ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << p3[i];
  }
  cout << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test064 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST064 tests PLANE_NORMAL2IMP_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double a;
  double b;
  double c;
  double d;
  double pn[DIM_NUM];
  double pp[DIM_NUM];

  cout << "\n";
  cout << "TEST064\n";
  cout << "  PLANE_NORMAL2IMP_3D puts a plane defined by\n";
  cout << "    point, normal form into implicit ABCD form.\n";
  cout << "\n";

  pp[0] = - 1.0;
  pp[1] =   0.0;
  pp[2] = - 1.0;

  pn[0] = - 0.2672612;
  pn[1] = - 0.5345225;
  pn[2] = - 0.8017837;

  cout << "\n";
  cout << "  Input:\n";
  cout << "\n";
  cout << "  PP = " << pp[0] << "  " << pp[1] << "  " << pp[2] << "\n";
  cout << "  PN = " << pn[0] << "  " << pn[1] << "  " << pn[2] << "\n";

  plane_normal2imp_3d ( pp, pn, &a, &b, &c, &d );

  cout << "\n";
  cout << "  Output:\n";
  cout << "\n";
  cout << "  (A,B,C,D) = " << a << "  " << b << "  " << c << "  " << d << "\n";

  pp[0] = - 16.0;
  pp[1] =    2.0;
  pp[2] =    4.0;

  pn[0] = - 0.2672612;
  pn[1] = - 0.5345225;
  pn[2] = - 0.8017837;

  plane_normal2imp_3d ( pp, pn, &a, &b, &c, &d );

  cout << "\n";
  cout << "  Input:\n";
  cout << "\n";
  cout << "  PP = " << pp[0] << "  " << pp[1] << "  " << pp[2] << "\n";
  cout << "  PN = " << pn[0] << "  " << pn[1] << "  " << pn[2] << "\n";
  cout << "\n";
  cout << "  Output:\n";
  cout << "\n";
  cout << "  (A,B,C,D) = " << a << "  " << b << "  " << c << "  " << d << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test065 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST065 tests POINTS_CENTROID_2D.
//
//  !....3&11....
//  !............
//  !............
//  X..9.........
//  !.....5......
//  !...........6
//  !.4.2...10...
//  !.....8...12.
//  V............
//  !..7.........
//  !......1.....
//  !............
//  !............
//  !----V----X--
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 12

  int centroid_index;
  double p[DIM_NUM*N] = {
     7.0,  3.0,
     4.0,  7.0,
     5.0, 13.0,
     2.0,  7.0,
     6.0,  9.0,
    12.0,  8.0,
     3.0,  4.0,
     6.0,  6.0,
     3.0, 10.0,
     8.0,  7.0,
     5.0, 13.0,
    10.0,  6.0 };

  cout << "\n";
  cout << "TEST065\n";
  cout << "  POINTS_CENTROID_2D computes the centroid of a\n";
  cout << "    discrete set of points.\n";

  r8mat_transpose_print ( DIM_NUM, N, p, "  The points:" );

  centroid_index = points_centroid_2d ( N, p );

  cout << "\n";
  cout << "  The centroid is point #: " << centroid_index << "\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test066 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST066 tests POINTS_COLIN_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double colin;
  int i;
  double *p1;
  double p1_test[DIM_NUM*TEST_NUM] = {
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0 };
  double *p2;
  double p2_test[DIM_NUM*TEST_NUM] = {
    10.0, 10.0,
     0.0, 1.0,
     1.0, 0.0 };
  double *p3;
  double p3_test[DIM_NUM*TEST_NUM] = {
      5.0, 4.99,
    100.0, 0.0,
      0.5, 0.86602539 };
  int test;

  cout << "\n";
  cout << "TEST066\n";
  cout << "  POINTS_COLIN_2D estimates the colinearity\n";
  cout << "    of three points.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p1_test + DIM_NUM * test;
    p2 = p2_test + DIM_NUM * test;
    p3 = p3_test + DIM_NUM * test;

    cout << "\n";

    if ( test == 1 )
    {
      cout << "  Points almost on a line: Expect tiny COLIN.\n";
    }
    else if ( test == 2 )
    {
      cout << "  Two points close, one far: Expect tiny COLIN.\n";
    }
    else if ( test == 3 )
    {
      cout << "  Equilateral triangle: Expect COLIN = 1.\n";
    }

    colin = points_colin_2d ( p1, p2, p3 );

    cout << "\n";
    cout << "  P1: ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << p1[i];
    }
    cout << "\n";
    cout << "  P2: ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << p2[i];
    }
    cout << "\n";
    cout << "  P3: ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << p3[i];
    }
    cout << "\n";

    cout << "\n";
    cout << "  Colinearity index = " << colin << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test068 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST068 tests the SPHERE_DISTANCE routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 6

  double dist1;
  double dist2;
  double dist3;
  string name[TEST_NUM] = {
    "Atlanta, Georgia  ",
    "North Pole        ",
    "South Pole        ",
    "Timbuktu          ",
    "San Antonio, Texas",
    "Savannah, Georgia " };
  int lat_d[TEST_NUM] =  { 33, 90, -90, 16, 29, 32 };
  int lat_m[TEST_NUM] =  { 11,  0,   0, 49, 25,  5 };
  int long_d[TEST_NUM] = { 82,  0,   0,  3, 98, 81 };
  int long_m[TEST_NUM] = { 34,  0,   0,  0, 30,  6 };
  double lat1;
  double lat2;
  double long1;
  double long2;
  double radius = 3957.0;
  int test1;
  int test2;

  cout << "\n";
  cout << "TEST068\n";
  cout << "  SPHERE_DISTANCE1, SPHERE_DISTANCE2 and SPHERE_DISTANCE3\n";
  cout << "  measure the distance between two points on a sphere.\n";
  cout << "\n";
  cout << "  All tests uses RADIUS = " << radius << "\n";
  cout << "  which is the radius of the earth in miles.\n";
  cout << "\n";

  for ( test1 = 0; test1 < TEST_NUM - 1; test1++ )
  {
    lat1 = dms_to_radians ( lat_d[test1], lat_m[test1], 0 );
    long1 = dms_to_radians ( long_d[test1], long_m[test1], 0 );

    cout << "\n";
    cout << "  Distance from " << name[test1] << "\n";

    for ( test2 = test1 + 1; test2 < TEST_NUM; test2++ )
    {
      lat2 = dms_to_radians ( lat_d[test2], lat_m[test2], 0 );
      long2 = dms_to_radians ( long_d[test2], long_m[test2], 0 );

      dist1 = sphere_distance1 ( lat1, long1, lat2, long2, radius );
      dist2 = sphere_distance2 ( lat1, long1, lat2, long2, radius );
      dist3 = sphere_distance3 ( lat1, long1, lat2, long2, radius );

      cout << "             to " << name[test2] 
           << "  " << dist1
           << "  " << dist2
           << "  " << dist3 << "\n";
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test0685 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0685 tests POLAR_TO_XY and XY_TO_POLAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 10

  double b;
  double c;
  double r;
  int seed;
  double t;
  int test;
  double xy1[2];
  double xy2[2];

  cout << "\n";
  cout << "TEST0685\n";
  cout << "  POLAR_TO_XY converts (R,Theta) to (X,Y);\n";
  cout << "  XY_TO_POLAR converts (X,Y) to (R,Theta).\n";
  cout << "\n";
  cout << "         X           Y     ===>  R           T   =>      X           Y\n";
  cout << "\n";

  b = -1.0;
  c = +1.0;
  seed = 123456789;

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    xy1[0] = r8_uniform ( b, c, seed );
    xy1[1] = r8_uniform ( b, c, seed );

    xy_to_polar ( xy1, &r, &t );
    polar_to_xy ( r, t, xy2 );

    cout << "  " << setw(10) << xy1[0]
         << "  " << setw(10) << xy1[1]
         << "  " << setw(10) << r
         << "  " << setw(10) << t
         << "  " << setw(10) << xy2[0]
         << "  " << setw(10) << xy2[1] << "\n";
  }

  return;
# undef TEST_NUM
}

//****************************************************************************80

void test0755 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0755 tests POLYGON_1_2D, POLYGON_X_2D, POLYGON_Y_2D, POLYGON_XX_2D etc.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 4

  double result;
  double v[DIM_NUM*N] = {
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0 };

  cout << "\n";
  cout << "TEST0755\n";
  cout << "  For a polygon in 2D:\n";
  cout << "  POLYGON_1_2D integrates 1\n";
  cout << "  POLYGON_X_2D integrates X\n";
  cout << "  POLYGON_Y_2D integrates Y\n";
  cout << "  POLYGON_XX_2D integrates X*X\n";
  cout << "  POLYGON_XY_2D integrates X*Y\n";
  cout << "  POLYGON_YY_2D integrates Y*Y\n";

  r8mat_transpose_print ( DIM_NUM, N, v, "  The polygon vertices:" );

  cout << "\n";
  cout << "  F(X,Y)    Integral\n";
  cout << "\n";

  result = polygon_1_2d ( N, v );
  cout << "    1  " << "  " << result << "\n";

  result = polygon_x_2d ( N, v );
  cout << "    X  " << "  " << result << "\n";

  result = polygon_y_2d ( N, v );
  cout << "    Y  " << "  " << result << "\n";

  result = polygon_xx_2d ( N, v );
  cout << "  X*X  " << "  " << result << "\n";

  result = polygon_xy_2d ( N, v );
  cout << "  X*Y  " << "  " << result << "\n";

  result = polygon_yy_2d ( N, v );
  cout << "  Y*Y  " << "  " << result << "\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test0757 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0757 tests POLYGON_ANGLES_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 6

  double *angle;
  int i;
  double v[DIM_NUM*N] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 1.0,
    3.0, 0.0,
    3.0, 2.0,
    1.0, 2.0 };

  cout << "\n";
  cout << "TEST0757\n";
  cout << "  For a polygon in 2D:\n";
  cout << "  POLYGON_ANGLES_2D computes the angles.\n";

  cout << "\n";
  cout << "  Number of polygonal vertices = " << N << "\n";

  r8mat_transpose_print ( DIM_NUM, N, v, "  The polygon vertices:" );

  angle = polygon_angles_2d ( N, v );

  cout << "\n";
  cout << "  Polygonal angles in degrees:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << radians_to_degrees ( angle[i] ) << "\n";
  }

  delete [] angle;

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test076 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST076 tests POLYGON_AREA_2D.
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
{
# define DIM_NUM 2
# define TEST_NUM 2

  double area;
  double area_exact;
  double area_exact_test[TEST_NUM] = { 2.0, 6.0 };
  int n;
  int n_test[TEST_NUM] = { 4, 8 };
  int test;
  double *v;
  double v1[DIM_NUM*4] = {
    1.0, 0.0,
    2.0, 1.0,
    1.0, 2.0,
    0.0, 1.0 };
  double v2[DIM_NUM*8] = {
    0.0, 0.0,
    3.0, 0.0,
    3.0, 3.0,
    2.0, 3.0,
    2.0, 1.0,
    1.0, 1.0,
    1.0, 2.0,
    0.0, 2.0 };

  cout << "\n";
  cout << "TEST076\n";
  cout << "  For a polygon in 2D:\n";
  cout << "  POLYGON_AREA_2D computes the area.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = n_test[test];
    area_exact = area_exact_test[test];
    v = new double[DIM_NUM*n];

    if ( test == 0 )
    {
      r8mat_copy ( DIM_NUM, n, v1, v );
    }
    else if ( test == 1 )
    {
      r8mat_copy ( DIM_NUM, n, v2, v );
    }

    cout << "\n";
    cout << "  Number of polygonal vertices = " << n << "\n";

    r8mat_transpose_print ( DIM_NUM, n, v, "  The polygon vertices:" );

    area = polygon_area_2d ( n, v );

    cout << "\n";
    cout << "  Exact area is        " << area_exact << "\n";
    cout << "  The computed area is " << area << "\n";

    delete [] v;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0765 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0765 tests POLYGON_AREA_2D_2.
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
{
# define DIM_NUM 2
# define TEST_NUM 2

  double area;
  double area_exact;
  double area_exact_test[TEST_NUM] = { 2.0, 6.0 };
  int n;
  int n_test[TEST_NUM] = { 4, 8 };
  int test;
  double *v;
  double v1[DIM_NUM*4] = {
    1.0, 0.0,
    2.0, 1.0,
    1.0, 2.0,
    0.0, 1.0 };
  double v2[DIM_NUM*8] = {
    0.0, 0.0,
    3.0, 0.0,
    3.0, 3.0,
    2.0, 3.0,
    2.0, 1.0,
    1.0, 1.0,
    1.0, 2.0,
    0.0, 2.0 };

  cout << "\n";
  cout << "TEST0765\n";
  cout << "  For a polygon in 2D:\n";
  cout << "  POLYGON_AREA_2D_2 computes the area.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = n_test[test];
    area_exact = area_exact_test[test];
    v = new double[DIM_NUM*n];

    if ( test == 0 )
    {
      r8mat_copy ( DIM_NUM, n, v1, v );
    }
    else if ( test == 1 )
    {
      r8mat_copy ( DIM_NUM, n, v2, v );
    }

    cout << "\n";
    cout << "  Number of polygonal vertices = " << n << "\n";

    r8mat_transpose_print ( DIM_NUM, n, v, "  The polygon vertices:" );

    area = polygon_area_2d_2 ( n, v );

    cout << "\n";
    cout << "  Exact area is        " << area_exact << "\n";
    cout << "  The computed area is " << area << "\n";

    delete [] v;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test078 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST078 tests POLYGON_AREA_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  double area;
  double area_exact;
  double area_exact_test[TEST_NUM] = { 2.4494898, 6.0 };
  int n;
  int n_test[TEST_NUM] = { 4, 8 };
  double normal[DIM_NUM];
  int test;
  double *v;
  double v1[DIM_NUM*4] = {
    1.0, 0.0, 0.0,
    2.0, 1.0, 1.0,
    1.0, 2.0, 1.0,
    0.0, 1.0, 0.0 };
  double v2[DIM_NUM*8] = {
    0.00000,  0.00000,  0.00000,
    2.62679,  1.26009, -0.715657,
    1.48153,  3.97300, -0.142512,
    0.605932, 3.55297,  0.0960401,
    1.36944,  1.74437, -0.286056,
    0.493842, 1.32433, -0.0475041,
    0.112090, 2.22864,  0.143544,
   -0.763505, 1.80861,  0.382097 };

  cout << "\n";
  cout << "TEST078\n";
  cout << "  For a polygon in 3D:\n";
  cout << "  POLYGON_AREA_3D computes the area;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    area_exact = area_exact_test[test];
    n = n_test[test];

    v = new double[DIM_NUM*n];

    if ( test == 0 )
    {
      r8mat_copy ( DIM_NUM, n, v1, v );
    }
    else if ( test == 1 )
    {
      r8mat_copy ( DIM_NUM, n, v2, v );
    }

    r8mat_transpose_print ( DIM_NUM, n, v, "  The polygon vertices:" );

    area = polygon_area_3d ( n, v, normal );

    cout << "\n";
    cout << "  Exact area is        " << area_exact << "\n";
    cout << "  The computed area is " << area << "\n";

    delete [] v;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0782 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0782 tests POLYGON_AREA_3D_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  double area;
  double area_exact;
  double area_exact_test[TEST_NUM] = { 2.4494898, 6.0 };
  int n;
  int n_test[TEST_NUM] = { 4, 8 };
  int test;
  double *v;
  double v1[DIM_NUM*4] = {
    1.0, 0.0, 0.0,
    2.0, 1.0, 1.0,
    1.0, 2.0, 1.0,
    0.0, 1.0, 0.0 };
  double v2[DIM_NUM*8] = {
    0.00000,  0.00000,  0.00000,
    2.62679,  1.26009, -0.715657,
    1.48153,  3.97300, -0.142512,
    0.605932, 3.55297,  0.0960401,
    1.36944,  1.74437, -0.286056,
    0.493842, 1.32433, -0.0475041,
    0.112090, 2.22864,  0.143544,
   -0.763505, 1.80861,  0.382097 };

  cout << "\n";
  cout << "TEST0782\n";
  cout << "  For a polygon in 3D:\n";
  cout << "  POLYGON_AREA_3D_2 computes the area;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    area_exact = area_exact_test[test];
    n = n_test[test];

    v = new double[DIM_NUM*n];

    if ( test == 0 )
    {
      r8mat_copy ( DIM_NUM, n, v1, v );
    }
    else if ( test == 1 )
    {
      r8mat_copy ( DIM_NUM, n, v2, v );
    }

    r8mat_transpose_print ( DIM_NUM, n, v, "  The polygon vertices:" );

    area = polygon_area_3d_2 ( n, v );

    cout << "\n";
    cout << "  Exact area is        " << area_exact << "\n";
    cout << "  The computed area is " << area << "\n";

    delete [] v;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test0784 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0784 tests POLYGON_CENTROID_2D and POLYGON_CENTROID_2D_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 4

  double *centroid;
  double v[DIM_NUM*N] = {
    1.0, 0.0,
    2.0, 1.0,
    1.0, 2.0,
    0.0, 1.0 };

  cout << "\n";
  cout << "TEST0784\n";
  cout << "  For a polygon in 2D:\n";
  cout << "  POLYGON_CENTROID_2D computes the centroid.\n";
  cout << "  POLYGON_CENTROID_2D_2 computes the centroid.\n";

  r8mat_transpose_print ( DIM_NUM, N, v, "  The polygon vertices:" );

  centroid = polygon_centroid_2d ( N, v );

  r8vec_print ( DIM_NUM, centroid, "  POLYGON_CENTROID_2D:" );

  delete [] centroid;

  centroid = polygon_centroid_2d_2 ( N, v );

  r8vec_print ( DIM_NUM, centroid, "  POLYGON_CENTROID_2D_2:" );

  delete [] centroid;

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test0786 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0786 tests POLYGON_CENTROID_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define N 4

  double *centroid;
  double v[DIM_NUM*N] = {
    1.0, 0.0, 0.0,
    2.0, 1.0, 1.0,
    1.0, 2.0, 1.0,
    0.0, 1.0, 0.0 };

  cout << "\n";
  cout << "TEST0786\n";
  cout << "  For a polygon in 3D:\n";
  cout << "  POLYGON_CENTROID_3D computes the centroid.\n";

  r8mat_transpose_print ( DIM_NUM, N, v, "  The polygon vertices:" );

  centroid = polygon_centroid_3d ( N, v );

  r8vec_print ( DIM_NUM, centroid, "  The centroid:" );

  delete [] centroid;

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test079 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST079 tests POLYGON_CONTAINS_POINT_2D and POLYGON_CONTAINS_POINT_2D_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 5
# define TEST_NUM 4

  bool inside1;
  bool inside2;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
    1.0,  1.0,
    3.0,  4.0,
    0.0,  2.0,
    0.5, -0.25 };
  int test;
  double v[DIM_NUM*N] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 1.0,
    1.0, 2.0,
    0.0, 2.0 };

  cout << "\n";
  cout << "TEST079\n";
  cout << "  POLYGON_CONTAINS_POINT_2D determines if \n";
  cout << "    a point is in a polygon.\n";
  cout << "  POLYGON_CONTAINS_POINT_2D_2 determines if\n";
  cout << "    a point is in a polygon.\n";

  r8mat_transpose_print ( DIM_NUM, N, v, "  The polygon vertices:" );

  cout << "\n";
  cout << "          P          In1  In2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + DIM_NUM * test;

    cout << "  " << setw(12) << p[0]
         << "  " << setw(12) << p[1];

    inside1 = polygon_contains_point_2d ( N, v, p );

    cout << "  " << setw(1) << inside1;

    inside2 = polygon_contains_point_2d_2 ( N, v, p );

    cout << "  " << setw(1) << inside2 << "\n";
  }

  return;
# undef DIM_NUM
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test080 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST080 tests POLYGON_DIAMETER_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 4

  double diameter;
  double diameter_exact = 2.0;
  double v[DIM_NUM*N] = {
    1.0, 0.0,
    2.0, 1.0,
    1.0, 2.0,
    0.0, 1.0 };

  cout << "\n";
  cout << "TEST080\n";
  cout << "  For a polygon in 2D:\n";
  cout << "  POLYGON_DIAMETER_2D computes the diameter;\n";

  r8mat_transpose_print ( DIM_NUM, N, v, "  The polygon vertices:" );

  diameter = polygon_diameter_2d ( N, v );

  cout << "\n";
  cout << "  Diameter ( computed ) " << diameter << "\n";
  cout << "  Diameter ( exact )    " << diameter_exact << "\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test0801 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0801 tests POLYGON_EXPAND_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 4

  double h;
  double v[DIM_NUM*N] = {
    1.0, 1.0,
    5.0, 1.0,
    2.0, 4.0,
    1.0, 3.0 };
  double *w;

  cout << "\n";
  cout << "TEST0801\n";
  cout << "  For a polygon in 2D:\n";
  cout << "  POLYGON_EXPAND_2D expands it by an amount H.\n";

  h = 0.5;

  r8mat_transpose_print ( DIM_NUM, N, v, "  The polygon vertices:" );

  cout << "\n";
  cout << "  The expansion amount H = " << h << "\n";

  w = polygon_expand_2d ( N, v, h );

  r8mat_transpose_print ( DIM_NUM, N, w, "  The expanded polygon:" );

  delete [] w;

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test0803 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0803 tests POLYGON_INRAD_DATA_2D, POLYGON_OUTRAD_DATA_2D, and POLYGON_SIDE_DATA_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  int n;
  double radin;
  double radout;
  double side;

  cout << "\n";
  cout << "TEST0803\n";
  cout << "  For a REGULAR polygon in 2D:\n";
  cout << "  the inradius, outradius and side are related.\n";
  cout << "  POLYGON_INRAD_DATA_2D uses the inradius;\n";
  cout << "  POLYGON_OUTRAD_DATA_2D uses the inradius;\n";
  cout << "  POLYGON_SIDE_DATA_2D uses the inradius;\n";

  for ( n = 3; n <= 5; n++ )
  {
    cout << "\n";
    cout << "  Number of polygonal sides = " << n << "\n";

    side = 1.0;

    cout << "\n";
    cout << "  Assuming SIDE = " << side << "\n";
    polygon_side_data_2d ( n, side, &area, &radin, &radout );
    cout << "    AREA =   " << area << "\n";
    cout << "    RADIN =  " << radin << "\n";
    cout << "    RADOUT = " << radout << "\n";

    cout << "\n";
    cout << "  Assuming RADIN = " << radin << "\n";
    polygon_inrad_data_2d ( n, radin, &area, &radout, &side );
    cout << "    AREA =   " << area << "\n";
    cout << "    RADOUT = " << radout << "\n";
    cout << "    SIDE =   " << side << "\n";

    cout << "\n";
    cout << "  Assuming RADOUT = " << radout << "\n";
    polygon_outrad_data_2d ( n, radout, &area, &radin, &side );
    cout << "    AREA =   " << area << "\n";
    cout << "    RADIN =  " << radin << "\n";
    cout << "    SIDE =   " << side << "\n";
  }
  return;
}
//****************************************************************************80

void test0805 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0805 tests POLYGON_DIAMETER_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 4

  double diameter;
  double diameter_exact = 2.0;
  double v[DIM_NUM*N] = {
    1.0, 0.0,
    2.0, 1.0,
    1.0, 2.0,
    0.0, 1.0 };

  cout << "\n";
  cout << "TEST0805\n";
  cout << "  For a polygon in 2D:\n";
  cout << "  POLYGON_DIAMETER_2D computes the diameter;\n";

  r8mat_transpose_print ( DIM_NUM, N, v, "  The polygonal vertices:" );

  diameter = polygon_diameter_2d ( N, v );

  cout << "\n";
  cout << "  Diameter ( computed ) " << diameter << "\n";
  cout << "  Diameter ( exact )    " << diameter_exact << "\n";

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void polygon_solid_angle_3d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_SOLID_ANGLE_3D_TEST tests POLYGON_SOLID_ANGLE_3D.
//
//  Discussion:
//
//    The polygonal vertices seemed to be stored by rows instead of
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double p[3];
  double solid_angle;
  int test;
  int test_num = 4;
  double *v;

  cout << "\n";
  cout << "POLYGON_SOLID_ANGLE_3D_TEST\n";
  cout << "  POLYGON_SOLID_ANGLE_3D computes the solid angle\n";
  cout << "  subtended by a planar polygon in 3D as viewed from\n";
  cout << "  a point P.\n";

  for ( test = 1; test <= test_num; test++ )
  {
//
//  One eighth of sphere surface, on the unit sphere surface.
//
    if ( test == 1 )
    {
      n = 3;

      v = new double[3*n];

      v[0+0*3] = 1.0;
      v[0+1*3] = 0.0;
      v[0+2*3] = 0.0;

      v[1+0*3] = 0.0;
      v[1+1*3] = 1.0;
      v[1+2*3] = 0.0;

      v[2+0*3] = 0.0;
      v[2+1*3] = 0.0;
      v[2+2*3] = 1.0;

      p[0] = 0.0;
      p[1] = 0.0;
      p[2] = 0.0;
    }
//
//  Reverse order of vertices.
//
    else if ( test == 2 )
    {
      n = 3;

      v = new double[3*n];

      v[0+0*3] = 1.0;
      v[0+1*3] = 0.0;
      v[0+2*3] = 0.0;

      v[1+0*3] = 0.0;
      v[1+1*3] = 0.0;
      v[1+2*3] = 1.0;

      v[2+0*3] = 0.0;
      v[2+1*3] = 1.0;
      v[2+2*3] = 0.0;

      p[0] = 0.0;
      p[1] = 0.0;
      p[2] = 0.0;
    }
//
//  One eighth of sphere surface, on the unit sphere surface, 
//  translated by (1,2,3).
//
    else if ( test == 3 )
    {
      n = 3;

      v = new double[3*n];

      v[0+0*3] = 2.0;
      v[0+1*3] = 1.0;
      v[0+2*3] = 1.0;

      v[1+0*3] = 2.0;
      v[1+1*3] = 3.0;
      v[1+2*3] = 2.0;

      v[2+0*3] = 3.0;
      v[2+1*3] = 3.0;
      v[2+2*3] = 4.0;

      p[0] = 1.0;
      p[1] = 2.0;
      p[2] = 3.0;
    }
//
//  One eighth of sphere surface, but on sphere of radius 2.
//
    else if ( test == 4 )
    {
      n = 3;

      v = new double[3*n];

      v[0+0*3] = 2.0;
      v[0+1*3] = 0.0;
      v[0+2*3] = 0.0;

      v[1+0*3] = 0.0;
      v[1+1*3] = 2.0;
      v[1+2*3] = 0.0;

      v[2+0*3] = 0.0;
      v[2+1*3] = 0.0;
      v[2+2*3] = 2.0;

      p[0] = 0.0;
      p[1] = 0.0;
      p[2] = 0.0;
    }

    cout << "\n";
    cout << "  TEST #" << test << "\n";
    cout << "\n";

    r8vec_print ( 3, p, "  The viewing point P:" );

    r8mat_transpose_print ( 3, n, v, "  The polygon vertices V:" );

    solid_angle = polygon_solid_angle_3d ( n, v, p );

    cout << "\n";
    cout << "  Solid angle subtended: " << solid_angle << "\n";

    delete [] v;
  }
  return;
}
//****************************************************************************80

void test081 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST081 tests POLYHEDRON_AREA_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define FACE_NUM 4
# define ORDER_MAX 3
# define NODE_NUM 4

  double area;
  double area_exact = 2.366025;
  double coord[DIM_NUM*NODE_NUM] =
  {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0 };
  int i;
  int j;
  int node[FACE_NUM * ORDER_MAX] =
  {
    2, 1, 0,
    0, 1, 3,
    0, 3, 2,
    1, 2, 3 };
  int order[FACE_NUM] = { 3, 3, 3, 3 };

  cout << "\n";
  cout << "TEST081\n";
  cout << "  For a polyhedron in 3D:\n";
  cout << "  POLYHEDRON_AREA_3D computes surface area;\n";
  cout << "\n";
  cout << "  Number of faces is " << FACE_NUM << "\n";

  i4vec_print ( FACE_NUM, order, "  Order of each face:" );

  cout << "\n";
  cout << "  Nodes per face:\n";
  cout << "\n";
  for ( i = 0; i < FACE_NUM; i++ )
  {
    cout << "  " << setw(4) << i;
    for ( j = 0; j < order[i]; j++ )
    {
      cout << "  " << setw(10) << node[i*ORDER_MAX+j];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Nodal coordinates:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( i = 0; i < 3; i++ )
    {
      cout << "  " << setw(10) << coord[i+j*3];
    }
    cout << "\n";
  }

  area = polyhedron_area_3d ( coord, ORDER_MAX, FACE_NUM, node, NODE_NUM,
    order );

  cout << "\n";
  cout << "  Surface area = " << area << "\n";
  cout << "  Exact area =   " << area_exact << "\n";

  return;
# undef DIM_NUM
# undef FACE_NUM
# undef ORDER_MAX
# undef NODE_NUM
}
//****************************************************************************80

void test082 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST082 tests POLYHEDRON_CENTROID_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define FACE_NUM 4
# define ORDER_MAX 3
# define NODE_NUM 4

  double *centroid;
  double centroid_exact[3] = { 0.25, 0.25, 0.25 };
  double coord[DIM_NUM*NODE_NUM] =
  {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0 };
  int i;
  int j;
  int node[FACE_NUM * ORDER_MAX] =
  {
    3, 2, 1,
    1, 2, 4,
    1, 4, 3,
    2, 3, 4 };
  int order[FACE_NUM] = { 3, 3, 3, 3 };

  cout << "\n";
  cout << "TEST082\n";
  cout << "  For a polyhedron in 3D:\n";
  cout << "  POLYHEDRON_CENTROID_3D computes the centroid;\n";
  cout << "\n";
  cout << "  Number of faces is " << FACE_NUM << "\n";

  i4vec_print ( FACE_NUM, order, "  Order of each face:" );

  cout << "\n";
  cout << "  Nodes per face:\n";
  cout << "\n";
  for ( j = 0; j < FACE_NUM; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( i = 0; i < order[j]; i++ )
    {
      cout << "  " << setw(10) << node[i+j*ORDER_MAX];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Nodal coordinates:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(10) << coord[i+j*DIM_NUM];
    }
    cout << "\n";
  }

  centroid = polyhedron_centroid_3d ( coord, ORDER_MAX, FACE_NUM, node,
    NODE_NUM, order );

  cout << "\n";
  cout << "  The computed centroid is "
    << setw(10) << centroid[0] << "  "
    << setw(10) << centroid[1] << "  "
    << setw(10) << centroid[2] << "\n";
  cout << "  The correct centroid is  "
    << setw(10) << centroid_exact[0] << "  "
    << setw(10) << centroid_exact[1] << "  "
    << setw(10) << centroid_exact[2] << "\n";

  delete [] centroid;

  return;
# undef DIM_NUM
# undef FACE_NUM
# undef ORDER_MAX
# undef NODE_NUM
}
//****************************************************************************80

void test0825 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0825 tests POLYHEDRON_CONTAINS_POINT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define FACE_NUM 4
# define FACE_ORDER_MAX 3
# define NODE_NUM 4
# define TEST_NUM 100

  double *c;
  int face_order[FACE_NUM] = { 3, 3, 3, 3};
  int face_point[FACE_ORDER_MAX*FACE_NUM] = {
    1, 2, 4,
    1, 3, 2,
    1, 4, 3,
    2, 3, 4 };
  bool inside1;
  bool inside2;
  double *p;
  int seed;
  int test;
  double v[DIM_NUM*NODE_NUM] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0 };

  cout << "\n";
  cout << "TEST0825\n";
  cout << "  POLYHEDRON_CONTAINS_POINT_3D determines if a point is \n";
  cout << "  inside a polyhedron.\n";
  cout << "\n";
  cout << "  We test this routine by using a tetrahedron as the polyhedron.\n";
  cout << "  For this shape, an independent check can be made,\n";
  cout << "  using barycentric coordinates.\n";
  cout << "\n";
  cout << "  We label these checks IN1 and IN2, and we expect them to agree.\n";

  r8mat_transpose_print ( DIM_NUM, NODE_NUM, v, "  The vertices:" );

  i4vec_print ( FACE_NUM, face_order, "  The face orders:" );

  i4mat_transpose_print ( FACE_ORDER_MAX, FACE_NUM, face_point,
    "  The nodes making each face:" );

  cout << "\n";
  cout << "      X           Y           Z      IN1 IN2\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    p = r8vec_uniform_01_new ( DIM_NUM, seed );

    inside1 = polyhedron_contains_point_3d ( NODE_NUM, FACE_NUM,
      FACE_ORDER_MAX, v, face_order, face_point, p );

    c = tetrahedron_barycentric_3d ( v, p );

    inside2 =  ( 0.0 <= c[0] ) && ( c[0] <= 1.0 ) &&
               ( 0.0 <= c[1] ) && ( c[1] <= 1.0 ) &&
               ( 0.0 <= c[2] ) && ( c[2] <= 1.0 ) &&
               ( 0.0 <= c[3] ) && ( c[3] <= 1.0 ) &&
               ( c[0] + c[1] + c[2] + c[3] <= 1.0 );

    cout << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << p[2]
         << "  " << setw(1) << inside1
         << "  " << setw(1) << inside2 << "\n";

    if ( inside1 != inside2 )
    {
      cout << "??? Disagreement!  Barycentric coordinates:\n";
      cout << "  " << setw(10) << c[0]
           << "  " << setw(10) << c[1]
           << "  " << setw(10) << c[2]
           << "  " << setw(10) << c[3] << "\n";
    }

    delete [] c;
    delete [] p;
  }
  return;
# undef DIM_NUM
# undef FACE_NUM
# undef FACE_ORDER_MAX
# undef NODE_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test083 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST083 tests POLYHEDRON_VOLUME_3D and POLYHEDRON_VOLUME_3D_2;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define FACE_NUM 4
# define NODE_NUM 4
# define ORDER_MAX 3

  double coord[DIM_NUM*NODE_NUM] =
  {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0 };
  int i;
  int j;
  int node[FACE_NUM * ORDER_MAX] =
  {
    2, 1, 0,
    0, 1, 3,
    0, 3, 2,
    1, 2, 3 };
  int order[FACE_NUM] = { 3, 3, 3, 3 };
  double volume_exact = 1.0 / 6.0;
  double volume1;
  double volume2;

  cout << "\n";
  cout << "TEST083\n";
  cout << "  For a polyhedron in 3D:\n";
  cout << "  POLYHEDRON_VOLUME_3D computes the volume;\n";
  cout << "  POLYHEDRON_VOLUME_3D_2 computes the volume;\n";
  cout << "\n";
  cout << "  Number of faces is " << FACE_NUM << "\n";

  i4vec_print ( FACE_NUM, order, "  Order of each face:" );

  cout << "\n";
  cout << "  Nodes per face:\n";
  cout << "\n";
  for ( j = 0; j < FACE_NUM; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( i = 0; i < order[j]; i++ )
    {
      cout << "  " << setw(10) << node[i+j*ORDER_MAX];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Nodal coordinates:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( i = 0; i < 3; i++ )
    {
      cout << "  " << setw(10) << coord[i+j*3];
    }
    cout << "\n";
  }

  volume1 = polyhedron_volume_3d ( coord, ORDER_MAX, FACE_NUM, node, NODE_NUM,
    order );

  volume2 = polyhedron_volume_3d_2 ( coord, ORDER_MAX, FACE_NUM, node, NODE_NUM,
    order );

  cout << "\n";
  cout << "  Volume (method 1) = " << volume1 << "\n";
  cout << "  Volume (method 2) = " << volume2 << "\n";
  cout << "  Volume (exact)    = " << volume_exact << "\n";

  return;
# undef DIM_NUM
# undef FACE_NUM
# undef NODE_NUM
# undef ORDER_MAX
}
//****************************************************************************80

void test084 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST084 tests POLYLINE_INDEX_POINT_ND and POLYLINE_ARCLENGTH_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 4

  int i;
  int j;
  double p[DIM_NUM*N] = {
    0.0, 0.0,
    1.0, 1.0,
    2.0, 0.0,
    0.0, 0.0 };
  double *pt;
  double *s;
  double t;

  t = 2.0;

  cout << "\n";
  cout << "TEST084\n";
  cout << "  POLYLINE_INDEX_POINT_ND finds a point on a \n";
  cout << "    polyline with given arclength.\n";
  cout << "  POLYLINE_ARCLENGTH_ND computes the arclength \n";
  cout << "    of the polyline, and its nodes.\n";
  cout << "\n";
  cout << "  The line we examine is defined by these points:\n";
//
//  The call to POLYLINE_ARCLENGTH_ND is just to help us believe the final result.
//
  s = polyline_arclength_nd ( DIM_NUM, N, p );

  cout << "\n";
  cout << "           P              Arclength(X,Y)\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p[i+j*DIM_NUM];
    }
    cout << "  " << setw(12) << s[j] << "\n";
  }
  cout << "\n";
  cout << "  We search for the point with coordinate " << t << "\n";

  pt = polyline_index_point_nd ( DIM_NUM, N, p, t );

  r8vec_print ( DIM_NUM, pt, "  The computed point:" );

  delete [] s;
  delete [] pt;

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test0844 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0844 tests POLYLINE_POINTS_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define NK 4
# define NT 13

  double pk[DIM_NUM*NK] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0,
    1.0, 2.0 };
  double *pt;

  cout << "\n";
  cout << "TEST0844\n";
  cout << "  POLYLINE_POINTS_ND computes points on a polyline.\n";

  r8mat_transpose_print ( DIM_NUM, NK, pk, "  The defining points:" );

  pt = polyline_points_nd ( DIM_NUM, NK, pk, NT );

  r8mat_transpose_print ( DIM_NUM, NT, pt, "  The computed points:" );

  delete [] pt;

  return;
# undef DIM_NUM
# undef NK
# undef NT
}
//****************************************************************************80

void test0845 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0845 tests POLYLOOP_ARCLENGTH_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 4

  int i;
  int j;
  int j2;
  double p[DIM_NUM*N] = {
    0.0, 0.0,
    1.0, 1.0,
    2.0, 0.0,
    0.0, 0.0 };
  double *s;

  cout << "\n";
  cout << "TEST0845\n";
  cout << "  POLYLOOP_ARCLENGTH_ND computes the arclength\n";
  cout << "    of the nodes of a polyloop.\n";

  s = polyloop_arclength_nd ( DIM_NUM, N, p );

  cout << "\n";
  cout << "           P            Arclength(P)\n";
  cout << "\n";

  for ( j = 0; j <= N; j++ )
  {
    j2 = i4_wrap ( j, 0, N-1 );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p[i+j2*DIM_NUM];
    }
    cout << "  " << s[j2] << "\n";
  }

  delete [] s;

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test0846 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0846 tests POLYLOOP_POINTS_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define NK 4
# define NT 12

  double pk[DIM_NUM*NK] = {
    0.0, 2.0,
    0.0, 0.0,
    1.0, 0.0,
    1.0, 2.0 };
  double *pt;

  cout << "\n";
  cout << "TEST0846\n";
  cout << "  POLYLOOP_POINTS_ND computes points on a polyloop.\n";

  r8mat_transpose_print ( DIM_NUM, NK, pk, "  The defining points:" );

  pt = polyloop_points_nd ( DIM_NUM, NK, pk, NT );

  r8mat_transpose_print ( DIM_NUM, NT, pt, "  The computed points:" );

  delete [] pt;

  return;
# undef DIM_NUM
# undef NK
# undef NT
}
//****************************************************************************80

void test085 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST085 tests PLANE_EXP_PRO3.
//
//  Discussion:
//
//    Projection is ( -1, 1, 1 ).
//    Projection is ( 4, 5, -8 ).
//    Projection is ( 0.33, 0.33, 0.33).
//    Projection is ( 5.33, -1.66, -2.66 ).
//    Projection is ( -1, 1, 1 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 5

  int i;
  double p1[DIM_NUM] = { 1.0, 0.0, 0.0 };
  double p2[DIM_NUM] = { 0.0, 1.0, 0.0 };
  double p3[DIM_NUM] = { 0.0, 0.0, 1.0 };
  double po[DIM_NUM*TEST_NUM] = {
    0.0,  2.0,  2.0,
    4.0,  5.0, -8.0,
    0.25, 0.25, 0.25,
    5.0, -2.0, -3.0,
   -2.0,  0.0,  0.0 };
  double pp[DIM_NUM*TEST_NUM];
  int test;

  cout << "\n";
  cout << "TEST085\n";
  cout << "  PLANE_EXP_PRO3 projects an object point\n";
  cout << "    orthographically into a plane.\n";
  cout << "\n";
  cout << "            PO                PP\n";
  cout << "\n";

  plane_exp_pro3 ( p1, p2, p3, TEST_NUM, po, pp );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(10) << po[i+test*DIM_NUM];
    }
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(10) << pp[i+test*DIM_NUM];
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test170 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST170 tests PROVEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 2

  double base[M*N] = {
    4.0, 3.0, 2.0, 1.0,
    1.0, 2.0, 3.0, 4.0 };
  double vecm[M] = { 1.0, 1.0, 1.0, 2.0 };
  double vecn[N];
  double vecnm[M];

  cout << "\n";
  cout << "TEST170\n";
  cout << "  PROVEC projects a vector onto a subspace.\n";

  r8mat_transpose_print ( M, N, base, "  Base vectors" );

  r8vec_print ( M, vecm, "  Vector to be projected:" );

  provec ( M, N, base, vecm, vecn, vecnm );

  r8vec_print ( N, vecn, "  Projected vector in BASE coordinates:" );

  r8vec_print ( M, vecnm, "  Projected vector in original coordinates:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test171 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST171 tests QUAD_AREA_2D and QUAD_AREA2_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2010
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double area;
  double quad[DIM_NUM*4] = {
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0 };

  cout << "\n";
  cout << "TEST171\n";
  cout << "  For a quadrilateral in 2D:\n";
  cout << "  QUAD_AREA_2D finds the area;\n";
  cout << "  QUAD_AREA2_2D finds the area;\n";

  r8mat_transpose_print ( DIM_NUM, 4, quad, "  The vertices:" );

  area = quad_area_2d ( quad );

  cout << "\n";
  cout << "  QUAD_AREA_2D area is  " << area << "\n";

  area = quad_area2_2d ( quad );

  cout << "  QUAD_AREA2_2D area is " << area << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test1712 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1712 tests QUAD_AREA_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double area1;
  double area2;
  int i;
  int j;
  double q[3*4] = {
    2.0, 2.0, 0.0, 
    0.0, 0.0, 0.0, 
    1.0, 1.0, 1.0, 
    3.0, 3.0, 1.0 };
  double t[3*3];

  cout << "\n";
  cout << "TEST1712\n";
  cout << "  For a quadrilateral in 3D:\n";
  cout << "  QUAD_AREA_3D finds the area.\n";

  r8mat_transpose_print ( 3, 4, q, "  The vertices:");

  area = quad_area_3d ( q );

  cout << "\n";
  cout << "  QUAD_AREA_3D area is     " << area << "\n";

  for ( j = 0; j < 3; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      t[i+j*3] = q[i+j*3];
    }
  }
  area1 = triangle_area_3d ( t );
  for ( j = 0; j < 2; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      t[i+j*3] = q[i+(j+2)*3];
    }
  }
  for ( i = 0; i < 3; i++ )
  {
    t[i+2*3] = q[i+0*3];
  }
  area2 = triangle_area_3d ( t );
  cout << "  Sum of TRIANGLE_AREA_3D: " << area1 + area2 << "\n";

  return;
}
//****************************************************************************80

void test1715 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1715 tests QUAD_CONTAINS_POINT_2D, QUAD_POINT_DIST_2D, QUAD_POINT_DIST_SIGNED_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 7

  double dist;
  double dist_signed;
  bool inside;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
     0.25,   0.25,
     0.75,   0.25,
     1.00,   1.00,
    11.00,   0.50,
     0.00,   0.50,
     0.50, -10.00,
     2.00,   2.00 };
  double q[DIM_NUM*4] = {
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0 };
  int test;

  cout << "\n";
  cout << "TEST1715\n";
  cout << "  For a quadrilateral in 2D:\n";
  cout << "  QUAD_AREA_2D finds the area;\n";
  cout << "  QUAD_CONTAINS_POINT_2D tells if a point is inside;\n";
  cout << "  QUAD_POINT_DIST_2D computes the distance.\n";
  cout << "  QUAD_POINT_DIST_SIGNED_2D computes signed distance.\n";

  r8mat_transpose_print ( DIM_NUM, 4, q, "  The vertices:" );

  cout << "\n";
  cout << "        P        Contains  Dist    Dist\n";
  cout << "                          Signed  Unsigned\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    inside = quad_contains_point_2d ( q, p );

    dist_signed = quad_point_dist_signed_2d ( q, p );

    dist = quad_point_dist_2d ( q, p );

    cout << "  " << setw(12) << p[0]
         << "  " << setw(12) << p[1]
         << "  " << setw(1) << inside
         << "  " << setw(12) << dist_signed
         << "  " << setw(12) << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test172 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST172 tests QUAT_CONJ, QUAT_INV, QUAT_MUL and QUAT_NORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 4

  int i;
  double q1[DIM_NUM] = { 2.0, 3.0, 4.0, 5.0 };
  double *q2;
  double *q3;

  cout << "\n";
  cout << "TEST172\n";
  cout << "  QUAT_CONJ conjugates a quaternion;\n";
  cout << "  QUAT_INV inverts a quaternion;\n";
  cout << "  QUAT_MUL multiplies quaternions.\n";
  cout << "  QUAT_NORM computes the norm.\n";
  cout << "\n";
  cout << "  Q1:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << q1[i];
  }
  cout << "\n";

  cout << "  Norm ( Q1 ) = " << quat_norm ( q1 ) << "\n";

  q2 = quat_conj ( q1 );

  cout << "  Q2 = conj ( Q1 ):";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << q2[i];
  }
  cout << "\n";

  q3 = quat_mul ( q1, q2 );

  cout << "  Q3 = Q1 * Q2:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << q3[i];
  }
  cout << "\n";

  delete [] q2;

  q2 = quat_inv ( q1 );

  cout << "  Q2 = inv ( Q1 ):";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << q2[i];
  }
  cout << "\n";

  delete [] q3;

  q3 = quat_mul ( q1, q2 );

  cout << "  Q3 = Q1 * Q2:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << q3[i];
  }
  cout << "\n";

  delete [] q2;
  delete [] q3;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test173 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST173 tests RADEC_DISTANCE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 6

  double dec1;
  double dec2;
  double *p1;
  double *p2;
  double p_test[DIM_NUM*TEST_NUM] = {
     1.0,  0.0,  0.0,
     0.0,  1.0,  0.0,
     0.0,  0.0,  1.0,
     1.0,  1.0,  1.0,
     5.0, -2.0, -1.0,
    -2.0, -2.0, -2.0 };
  double ra1;
  double ra2;
  int test1;
  int test2;
  double theta;
  double theta_deg;

  cout << "\n";
  cout << "TEST173\n";
  cout << "  RADEC_DISTANCE_3D computes the angular separation\n";
  cout << "    between two points on a sphere described in terms of\n";
  cout << "    right ascension and declination.\n";
  cout << "\n";
  cout << "     RA1       DEC1      RA2       DEC2    Radians  Degrees\n";
  cout << "\n";

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    p1 = p_test + test1 * DIM_NUM;

    xyz_to_radec ( p1, &ra1, &dec1 );

    for ( test2 = test1 + 1; test2 < TEST_NUM; test2++ )
    {
      p2 = p_test + test2 * DIM_NUM;

      xyz_to_radec ( p2, &ra2, &dec2 );
      theta = radec_distance_3d ( ra1, dec1, ra2, dec2 );
      theta_deg = radians_to_degrees ( theta );

      cout << "  " << setw(8) << ra1
           << "  " << setw(8) << dec1
           << "  " << setw(8) << ra2
           << "  " << setw(8) << dec2
           << "  " << setw(8) << theta
           << "  " << setw(8) << theta_deg << "\n";
    }
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test174 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST174 tests RADEC_TO_XYZ and XYZ_TO_RADEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 6

  double dec;
  int i;
  double *p1;
  double *p2;
  double p_test[DIM_NUM*TEST_NUM] = {
     1.0,  0.0,  0.0,
     0.0,  1.0,  0.0,
     0.0,  0.0,  1.0,
     1.0,  1.0,  1.0,
     5.0, -2.0, -1.0,
    -2.0, -2.0, -2.0 };
  double ra;
  int test;

  cout << "\n";
  cout << "TEST174\n";
  cout << "  RADEC_TO_XYZ converts XYZ to RADEC coordinates.\n";
  cout << "  XYZ_TO_RADEC converts RADEC to XYZ coordinates.\n";
  cout << "\n";
  cout << "          P1          RA     DEC           P2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p_test + test * DIM_NUM;

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(7) << p1[i];
    }

    xyz_to_radec ( p1, &ra, &dec );

    cout << setw(7) << ra;
    cout << setw(7) << dec;

    p2 = radec_to_xyz ( ra, dec );

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(7) << p2[i];
    }
    cout << "\n";

    delete [] p2;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test1745 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1745 tests R8MAT_SOLVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define RHS_NUM 2

  double a[N*(N+RHS_NUM)] = {
     1.0,  4.0,  7.0,
     2.0,  5.0,  8.0,
     3.0,  6.0,  0.0,
    14.0, 32.0, 23.0,
     7.0, 16.0,  7.0 };
  int i;
  int info;
  int j;

  cout << "\n";
  cout << "TEST1745\n";
  cout << "  R8MAT_SOLVE solves linear systems.\n";
//
//  Print out the matrix to be inverted.
//
  r8mat_print ( N, N+RHS_NUM, a, "  The linear system:" );
//
//  Solve the systems.
//
  info = r8mat_solve ( N, RHS_NUM, a );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  The input matrix was singular.\n";
    cout << "  The solutions could not be computed.\n";
    return;
  }

  cout << "\n";
  cout << "  The computed solutions:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    for ( j = N; j < N+RHS_NUM; j++ )
    {
      cout << "  " << setw(10) << a[i+j*N];
    }
    cout << "\n";
  }

  return;
# undef N
# undef RHS_NUM
}
//****************************************************************************80********

void test1746 ( )

//****************************************************************************80********
//
//  Purpose:
//
//    TEST1746 tests R8MAT_INVERSE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a[3*3] = {
    3.0, 2.0, 0.0,
    2.0, 2.0, 1.0,
    1.0, 1.0, 1.0 };
  double *b;
  int i;

  cout << "\n";
  cout << "TEST1746\n";
  cout << "  R8MAT_INVERSE_3D inverts a 3 by 3 matrix.\n";

  cout << "\n";
  cout << "  Matrix A:\n";
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << setw(10) << a[i+0*3]
         << "  " << setw(10) << a[i+1*3]
         << "  " << setw(10) << a[i+2*3] << "\n";
  }

  b = r8mat_inverse_3d ( a );

  cout << "\n";
  cout << "  Inverse matrix B:\n";
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << setw(10) << b[i+0*3]
         << "  " << setw(10) << b[i+1*3]
         << "  " << setw(10) << b[i+2*3] << "\n";
  }

  delete [] b;

  return;
}
//****************************************************************************80

void test175 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST175 tests ROTATION_AXIS_VECTOR_3D, ROTATION_MAT_VECTOR_3D, ROTATION_QUAT_VECTOR_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double a[DIM_NUM*DIM_NUM];
  double angle = 1.159804;
  double axis[DIM_NUM] = { 0.2361737, -0.8814124, -0.4090649 };
  double q[4];
  double v[DIM_NUM] = { 1.0, 4.0, 10.0 };
  double w[DIM_NUM];

  cout << "\n";
  cout << "TEST175\n";
  cout << "  ROTATION_AXIS_VECTOR_3D applies an axis\n";
  cout << "    rotation to a vector;\n";
  cout << "  ROTATION_MAT_VECTOR_3D applies a matrix\n";
  cout << "    rotation to a vector.\n";
  cout << "  ROTATION_QUAT_VECTOR_3D applies a quaternion\n";
  cout << "    rotation to a vector.\n";

  r8vec_print ( DIM_NUM, v, "  The vector:" );

  r8vec_print ( DIM_NUM, axis, "  The rotation axis:" );

  cout << "\n";
  cout << "  The rotation angle is " << angle << "\n";

  rotation_axis_vector_3d ( axis, angle, v, w );

  r8vec_print ( DIM_NUM, w, "  The rotated vector:" );

  rotation_axis2mat_3d ( axis, angle, a );

  r8mat_print ( DIM_NUM, DIM_NUM, a, "  The rotation matrix:" );

  rotation_mat_vector_3d ( a, v, w );

  r8vec_print ( DIM_NUM, w, "  The rotated vector:" );

  rotation_axis2quat_3d ( axis, angle, q );

  r8vec_print ( 4, q, "  The rotation quaternion:" );

  rotation_quat_vector_3d ( q, v, w );

  r8vec_print ( DIM_NUM, w, "  The rotated vector:" );
//
//  Another test of ROTATION_AXIS_VECTOR_3D with an axis vector
//  that does not have unit length.
//
  v[0] = 1.0;
  v[1] = 1.0;
  v[2] = 1.0;

  r8vec_print ( DIM_NUM, v, "  The vector:" );

  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 2.0;

  r8vec_print ( DIM_NUM, axis, "  The rotation axis:" );

  angle = 90.0;
  angle = degrees_to_radians ( angle );

  cout << "\n";
  cout << "  The rotation angle is " << angle << "\n";

  rotation_axis_vector_3d ( axis, angle, v, w );

  r8vec_print ( DIM_NUM, w, "  The rotated vector:" );

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test176 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST176 tests ROTATION_AXIS2MAT_3D and ROTATION_MAT2AXIS_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double a[DIM_NUM*DIM_NUM] = {
    0.43301269, -0.5,        0.75,
    0.25,        0.86602539, 0.43301269,
   -0.86602539,  0.0,        0.5 };
  double b[DIM_NUM*DIM_NUM];
  double angle;
  double axis[DIM_NUM];

  cout << "\n";
  cout << "TEST176\n";
  cout << "  ROTATION_MAT2AXIS_3D computes a rotation axis\n";
  cout << "    and angle from a rotation matrix.\n";
  cout << "  ROTATION_AXIS2MAT_3D computes a rotation matrix\n";
  cout << "    from a rotation axis and angle.\n";

  r8mat_print ( DIM_NUM, DIM_NUM, a, "  The rotation matrix:" );

  rotation_mat2axis_3d ( a, axis, &angle );

  r8vec_print ( DIM_NUM, axis, "  The rotation axis:" );

  cout << "\n";
  cout << "  The rotation angle is " << angle << "\n";

  rotation_axis2mat_3d ( axis, angle, b );

  r8mat_print ( DIM_NUM, DIM_NUM, b, "  The rotation matrix:" );

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test177 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST177 tests ROTATION_AXIS2QUAT_3D, ROTATION_QUAT2AXIS_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double angle = 1.159804;
  double axis[DIM_NUM] = { 0.2361737, -0.8814124, -0.4090649 };
  double q[4];

  cout << "\n";
  cout << "TEST177\n";
  cout << "  ROTATION_QUAT2AXIS_3D computes a rotation axis\n";
  cout << "    and angle from a rotation quaternion.\n";
  cout << "  ROTATION_AXIS2QUAT_3D computes a rotation\n";
  cout << "    quaternion from a rotation axis and angle.\n";

  r8vec_print ( DIM_NUM, axis, "  Rotation axis:" );

  cout << "  Rotation angle is " << angle << "\n";

  rotation_axis2quat_3d ( axis, angle, q );

  r8vec_print ( 4, q, "  The rotation quaternion:" );

  rotation_quat2axis_3d ( q, axis, &angle );

  r8vec_print ( DIM_NUM, axis, "  The rotation axis:" );

  cout << "\n";
  cout << "  The rotation angle is " << angle << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test178 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST178 tests ROTATION_MAT2QUAT_3D and ROTATION_QUAT2MAT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double a[DIM_NUM*DIM_NUM] = {
    0.43301269, -0.5,        0.75,
    0.25,        0.86602539, 0.43301269,
   -0.86602539,  0.0,        0.5 };
  double b[DIM_NUM*DIM_NUM];
  double q[4];

  cout << "\n";
  cout << "TEST178\n";
  cout << "  ROTATION_MAT2QUAT_3D computes a rotation\n";
  cout << "    quaternion from a rotation matrix.\n";
  cout << "  ROTATION_QUAT2MAT_3D computes a rotation matrix\n";
  cout << "    from a rotation quaternion.\n";

  r8mat_print ( DIM_NUM, DIM_NUM, a, "  The rotation matrix:" );

  rotation_mat2quat_3d ( a, q );

  r8vec_print ( 4, q, "  The rotation quaternion:" );

  rotation_quat2mat_3d ( q, b );

  r8mat_print ( DIM_NUM, DIM_NUM, b, "  The rotation matrix:" );

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test1787 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1787 tests DGE_FA and DGE_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double alu[N*N];
  double b[N];
  int i;
  int info;
  int j;
  int job;
  int pivot[N];
  int seed;
  double x[N];

  cout << "\n";
  cout << "TEST1787\n";
  cout << "  DGE_FA factors a general linear system,\n";
  cout << "  DGE_SL solves a factored system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  seed = 123456789;
  a = r8mat_uniform_01_new ( N, N, seed );
  r8mat_print ( N, N, a, "  Matrix A:" );
//
//  Set the desired solution.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  r8vec_print ( N, x, "  Desired solution vector:" );
//
//  Compute the corresponding right hand side.
//
  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a[i+j*N] * x[j];
    }
  }
  r8vec_print ( N, b, "  Right hand side vector:" );
//
//  Make a copy of the matrix.
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      alu[i+j*N] = a[i+j*N];
    }
  }
//
//  Factor the matrix.
//
  info = dge_fa ( N, alu, pivot );

  r8mat_print ( N, N, alu, "  Factored matrix ALU:" );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Fatal error!\n";
    cout << "  DGE_FA declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  job = 0;
  dge_sl ( N, alu, pivot, b, job );

  r8vec_print ( N, b, "  Solution: (Should be 1, 2, 3,...)" );
//
//  Set another the desired solution.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = 1.0;
  }
//
//  Compute the corresponding right hand side.
//
  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a[i+j*N] * x[j];
    }
  }
//
//  Solve the system
//
  job = 0;
  dge_sl ( N, alu, pivot, b, job );

  r8vec_print ( N, b, "  Solution: (Should be 1, 1, 1,...)" );
//
//  Set the desired solution to a problem involving the transposed matrix.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
//
//  Compute the corresponding right hand side using the transposed matrix.
//
  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a[j+i*N] * x[j];
    }
  }
//
//  Solve the transposed system.
//
  job = 1;
  dge_sl ( N, alu, pivot, b, job );

  r8vec_print ( N, b,
    "  Solution of transposed system: (Should be 1, 2, 3,...)" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test1893 ( )

//****************************************************************************80
//
//  Purpose:
//
//   TEST1893 tests RTP_TO_XYZ and XYZ_TO_RTP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    221 July 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a = -2.0;
  double b =  3.0;
  double phi;
  double r;
  int seed;
  int test;
  int test_num = 5;
  double theta;
  double *xyz1;
  double xyz2[3];

  cout << "\n";
  cout << "TEST1893\n";
  cout << "  RTP_TO_XYZ converts XYZ to (R,Theta,Phi) coordinates.\n";
  cout << "  XYZ_TO_RTP converts (R,Theta,Phi) to XYZ coordinates.\n";
  cout << "\n";
  cout << "      X1     Y1     Z1     R    THETA    PHI    X2     Y2     Z2\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    xyz1 = r8vec_uniform_new ( 3, a, b, seed );

    xyz_to_rtp ( xyz1, &r, &theta, &phi );
    rtp_to_xyz ( r, theta, phi, xyz2 );

    cout << "  " << setw(7) << xyz1[0]
         << "  " << setw(7) << xyz1[1]
         << "  " << setw(7) << xyz1[1]
         << "  " << setw(7) << r
         << "  " << setw(7) << theta
         << "  " << setw(7) << phi
         << "  " << setw(7) << xyz2[0]
         << "  " << setw(7) << xyz2[1]
         << "  " << setw(7) << xyz2[2] << "\n";

    delete [] xyz1;
  }

  return;
}
//****************************************************************************80

void test036 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST036 tests SEGMENT_CONTAINS_POINT_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  double p;
  double p_test[TEST_NUM] = { 3.0,   7.5, 20.0,  5.0 };
  double p1;
  double p1_test[TEST_NUM] = { 2.0,  10.0,  8.0, 88.0 };
  double p2;
  double p2_test[TEST_NUM] = { 6.0, -10.0, 10.0, 88.0 };
  double t;
  int test;

  cout << "\n";
  cout << "TEST036\n";
  cout << "  SEGMENT_CONTAINS_POINT_1D determines if a point\n";
  cout << "    lies within a line segment in 1D.\n";
  cout << "\n";
  cout << "       P1     P       T\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p1_test[test];
    p2 = p2_test[test];
    p = p_test[test];

    segment_contains_point_1d ( p1, p2, p, &t );

    cout << "  " << setw(7) << p1
         << "  " << setw(7) << p2
         << "  " << setw(7) << p
         << "  " << setw(12) << t << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test0365 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0365 tests SEGMENT_POINT_DIST_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double dist;
  int i;
  double *p;
  double *p1;
  double *p2;
  int seed = 123456789;
  int test;
  int test_num = 3;

  cout << "\n";
  cout << "TEST0365\n";
  cout << "  SEGMENT_POINT_DIST_2D computes the distance\n";
  cout << "    between a line segment and point in 2D\n";

  for ( test = 1; test <= test_num; test++ )
  {
    p1 = r8vec_uniform_01_new ( DIM_NUM, seed );
    p2 = r8vec_uniform_01_new ( DIM_NUM, seed );
    p = r8vec_uniform_01_new ( DIM_NUM, seed );

    dist = segment_point_dist_2d ( p1, p2, p );

    cout << "\n";
    cout << "  TEST = " << setw(2) << test << "\n";
    cout << "  P1 =   ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p1[i];
    }
    cout << "\n";
    cout << "  P2 =   ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p2[i];
    }
    cout << "\n";
    cout << "  P =    ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p[i];
    }
    cout << "\n";
    cout << "  DIST = " << setw(12) << dist << "\n";

    delete [] p1;
    delete [] p2;
    delete [] p;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0366 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0366 tests SEGMENT_POINT_DIST_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double dist;
  int i;
  double *p;
  double *p1;
  double *p2;
  int seed = 123456789;
  int test;
  int test_num = 3;

  cout << "\n";
  cout << "TEST0366\n";
  cout << "  SEGMENT_POINT_DIST_3D computes the distance\n";
  cout << "    between a line segment and point in 3D\n";

  for ( test = 1; test <= test_num; test++ )
  {
    p1 = r8vec_uniform_01_new ( DIM_NUM, seed );
    p2 = r8vec_uniform_01_new ( DIM_NUM, seed );
    p = r8vec_uniform_01_new ( DIM_NUM, seed );

    dist = segment_point_dist_3d ( p1, p2, p );

    cout << "\n";
    cout << "  TEST = " << setw(2) << test << "\n";
    cout << "  P1 =   ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p1[i];
    }
    cout << "\n";
    cout << "  P2 =   ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p2[i];
    }
    cout << "\n";
    cout << "  P =    ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p[i];
    }
    cout << "\n";
    cout << "  DIST = " << setw(12) << dist << "\n";

    delete [] p1;
    delete [] p2;
    delete [] p;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0367 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0367 tests SEGMENT_POINT_NEAR_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double dist;
  int i;
  double *p;
  double *p1;
  double *p2;
  double pn[DIM_NUM];
  int seed = 123456789;
  double t;
  int test;
  int test_num = 3;

  cout << "\n";
  cout << "TEST0367\n";
  cout << "  SEGMENT_POINT_NEAR_2D computes the nearest point\n";
  cout << "    on a line segment to a point in 2D\n";

  for ( test = 1; test <= test_num; test++ )
  {
    p1 = r8vec_uniform_01_new ( DIM_NUM, seed );
    p2 = r8vec_uniform_01_new ( DIM_NUM, seed );
    p = r8vec_uniform_01_new ( DIM_NUM, seed );

    segment_point_near_2d ( p1, p2, p, pn, &dist, &t );

    cout << "\n";
    cout << "  TEST = " << setw(2) << test << "\n";
    cout << "  P1 =   ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p1[i];
    }
    cout << "\n";
    cout << "  P2 =   ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p2[i];
    }
    cout << "\n";
    cout << "  P =    ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p[i];
    }
    cout << "\n";
    cout << "  PN =   ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << pn[i];
    }
    cout << "\n";
    cout << "  DIST = " << setw(12) << dist << "\n";
    cout << "  T =    " << setw(12) << t << "\n";
    delete [] p1;
    delete [] p2;
    delete [] p;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test0368 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0368 tests SEGMENT_POINT_NEAR_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double dist;
  int i;
  double *p;
  double *p1;
  double *p2;
  double pn[DIM_NUM];
  int seed = 123456789;
  double t;
  int test;
  int test_num = 3;

  cout << "\n";
  cout << "TEST0368\n";
  cout << "  SEGMENT_POINT_NEAR_3D computes the nearest point\n";
  cout << "    on a line segment to a point in 3D\n";

  for ( test = 1; test <= test_num; test++ )
  {
    p1 = r8vec_uniform_01_new ( DIM_NUM, seed );
    p2 = r8vec_uniform_01_new ( DIM_NUM, seed );
    p = r8vec_uniform_01_new ( DIM_NUM, seed );

    segment_point_near_3d ( p1, p2, p, pn, &dist, &t );

    cout << "\n";
    cout << "  TEST = " << setw(2) << test << "\n";
    cout << "  P1 =   ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p1[i];
    }
    cout << "\n";
    cout << "  P2 =   ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p2[i];
    }
    cout << "\n";
    cout << "  P =    ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << p[i];
    }
    cout << "\n";
    cout << "  PN =   ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << pn[i];
    }
    cout << "\n";
    cout << "  DIST = " << setw(12) << dist << "\n";
    cout << "  T    = " << setw(12) << t << "\n";
    delete [] p1;
    delete [] p2;
    delete [] p;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test037 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST037 tests SEGMENT_POINT_NEAR_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double dist;
  double p[DIM_NUM];
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double pn[DIM_NUM];
  double t;

  cout << "\n";
  cout << "TEST037\n";
  cout << "  SEGMENT_POINT_NEAR_3D computes the nearest point on a\n";
  cout << "  line segment, to a given point, in 3 space.\n";
  cout << "\n";
  cout << "  Case  T  Distance    PN\n";
  cout << "\n";
//
//  Case 1, point is nearest end of segment.
//
//  LS: (2,3,0) + t * (2,1,0) for t = 0 to 3.
//  P (11,6,4)
//  Distance is 5.
//
  p1[0] = 2.0;
  p1[1] = 3.0;
  p1[2] = 0.0;

  p2[0] = 8.0;
  p2[1] = 6.0;
  p2[2] = 0.0;

  p[0] = 11.0;
  p[1] =  6.0;
  p[2] =  4.0;

  segment_point_near_3d ( p1, p2, p, pn, &dist, &t );

  cout << "  " << setw(6)  << 1
       << "  " << setw(10) << t
       << "  " << setw(10) << dist
       << "  " << setw(10) << pn[0]
       << "  " << setw(10) << pn[1]
       << "  " << setw(10) << pn[2] << "\n";
//
//  Case 2, point is nearest interior point of segment.
//
//  LS: (2,3,0) + t * (2,1,0) for t = 0 to 3.
//  P (4,4,1)
//  Distance is 1.
//
  p1[0] = 2.0;
  p1[1] = 3.0;
  p1[2] = 0.0;

  p2[0] = 8.0;
  p2[1] = 6.0;
  p2[2] = 0.0;

  p[0] =  4.0;
  p[1] =  4.0;
  p[2] =  1.0;

  segment_point_near_3d ( p1, p2, p, pn, &dist, &t );

  cout << "  " << setw(6)  << 2
       << "  " << setw(10) << t
       << "  " << setw(10) << dist
       << "  " << setw(10) << pn[0]
       << "  " << setw(10) << pn[1]
       << "  " << setw(10) << pn[2] << "\n";
//
//  Case 3, point is on the line.
//
//  LS: (2,3,0) + t * (2,1,0) for t = 0 to 3.
//  P (6,5,0)
//  Distance is 0.
//
  p1[0] = 2.0;
  p1[1] = 3.0;
  p1[2] = 0.0;

  p2[0] = 8.0;
  p2[1] = 6.0;
  p2[2] = 0.0;

  p[0] =  6.0;
  p[1] =  5.0;
  p[2] =  0.0;

  segment_point_near_3d ( p1, p2, p, pn, &dist, &t );

  cout << "  " << setw(6)  << 3
       << "  " << setw(10) << t
       << "  " << setw(10) << dist
       << "  " << setw(10) << pn[0]
       << "  " << setw(10) << pn[1]
       << "  " << setw(10) << pn[2] << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test1788 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17888 tests SIMPLEX_LATTICE_LAYER_POINT_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  int *c;
  int i;
  int j;
  int layer;
  bool more;
  int n;
  int n_test[TEST_NUM] = { 1, 2, 3, 4 };
  int test;
  int *v;

  cout << "\n";
  cout << "TEST1788\n";
  cout << "  SIMPLEX_LATTICE_LAYER_POINT_NEXT returns the next\n";
  cout << "  point in an N-dimensional simplex lattice layer defined by:\n";
  cout << "\n";
  cout << "    C(N+1) - 1 <= X[0]/C[0] + X(2)/C(2) + ... + X(N)/C(N) <= C(N+1).\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = n_test[test];
    c = new int[n+1];
    v = new int[n];

    for ( i = 0; i < n + 1; i++ )
    {
      c[i] = i + 2;
    }
    for ( i = 0; i < n; i++ )
    {
      v[i] = 0;
    }

    cout << "\n";
    cout << "  N = " << n << "\n";
    cout << "  C = ";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(4) << c[i];
    }
    cout << "\n";
    cout << "\n";

    for ( layer = 0; layer <= 2; layer++ )
    {
      cout << "\n";
      cout << "  Layer " << layer << "\n";
      cout << "\n";

      c[n] = layer;
      more = false;
      i = 0;

      for ( ; ; )
      {
        simplex_lattice_layer_point_next ( n, c, v, &more );
        if ( !more )
        {
          cout << "  No more.\n";
          break;
        }
        i = i + 1;
        cout << "  " << setw(4) << i;
        for ( j = 0; j < n; j++ )
        {
          cout << "  " << setw(4) << v[j];
        }
        cout << "\n";
      }
    }
    delete [] c;
    delete [] v;
  }

  return;
# undef TEST_NUM
}

//****************************************************************************80

void test1789 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1789 tests SIMPLEX_LATTICE_POINT_NEXT.
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
{
# define TEST_NUM 4

  int *c;
  int i;
  int j;
  bool more;
  int n;
  int n_test[TEST_NUM] = { 1, 2, 3, 4 };
  int test;
  int *v;

  cout << "\n";
  cout << "TEST1789\n";
  cout << "  SIMPLEX_LATTICE_POINT_NEXT returns the next lattice\n";
  cout << "  point in an N-dimensional simplex defined by:\n";
  cout << "\n";
  cout << "    0 <= X(1)/C(1) + X(2)/C(2) + ... + X(N)/C(N) <= C(N+1).\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = n_test[test];

    c = new int[n+1];
    v = new int[n];

    for ( i = 0; i < n + 1; i++ )
    {
      c[i] = n + 1 - i;
    }
    for ( i = 0; i < n; i++ )
    {
      v[i] = 0;
    }
    more = false;

    cout << "\n";
    cout << "  N = " << n << "\n";
    cout << "  C = ";
    for ( i = 0; i < n + 1; i++ )
    {
      cout << "  " << setw(4) << c[i];
    }
    cout << "\n";
    cout << "\n";

    i = 0;

    for ( ; ; )
    {
      simplex_lattice_point_next ( n, c, v, &more );
      if ( !more )
      {
        cout << "  No more.\n";
        break;
      }
      i = i + 1;
      cout << "  " << setw(4) << i;
      for ( j = 0; j < n; j++ )
      {
        cout << "  " << setw(4) << v[j];
      }
      cout << "\n";
    }
    delete [] c;
    delete [] v;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test179 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST179 tests SOCCER_SIZE_3D, SOCCER_SHAPE_3D, SHAPE_PRINT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int edge_num;
  int face_num;
  int *face_order;
  int face_order_max;
  int *face_point;
  int point_num;
  double *point_coord;

  cout << "\n";
  cout << "TEST179\n";
  cout << "  For the truncated icosahedron, or soccer ball,\n";
  cout << "  SOCCER_SIZE_3D returns dimension information;\n";
  cout << "  SOCCER_SHAPE_3D returns face and order information.\n";
  cout << "  SHAPE_PRINT_3D prints this information.\n";
//
//  Get the data sizes.
//
  soccer_size_3d ( &point_num, &edge_num, &face_num, &face_order_max );

  cout << "\n";
  cout << "    Number of vertices: " << point_num << "\n";
  cout << "    Number of edges   : " << edge_num << "\n";
  cout << "    Number of faces   : " << face_num << "\n";
  cout << "    Maximum face order: " << face_order_max << "\n";
//
//  Make room for the data.
//
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];
  point_coord = new double[DIM_NUM*point_num];
//
//  Get the data.
//
  soccer_shape_3d ( point_num, face_num, face_order_max, point_coord,
    face_order, face_point );
//
//  Print the data.
//
  shape_print_3d ( point_num, face_num, face_order_max,
    point_coord, face_order, face_point );

  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test180 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST180 tests SORT_HEAP_EXTERNAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int a[N];
  int i;
  int indx;
  int isgn;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST180\n";
  cout << "  SORT_HEAP_EXTERNAL sorts objects externally.\n";

  indx = 0;
  i = 0;
  j = 0;
  isgn = 0;
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i] = i4_uniform ( 1, N, seed );
  }

  i4vec_print ( N, a, "  Unsorted array:" );
//
//void the sort routine over and over.
//
  for ( ;; )
  {
    sort_heap_external ( N, &indx, &i, &j, isgn );
//
//  If the return value of INDX is negative, we're asked to compare
//  array elements I and J;
//
    if ( indx < 0 )
    {
      if ( a[i-1] <= a[j-1] )
      {
        isgn = -1;
      }
      else
      {
        isgn = 1;
      }
    }
//
//  ...and if the return value of INDX is positive, we're asked to switch
//  array elements I and J;
//
    else if ( 0 < indx )
    {
      i4_swap ( &a[i-1], &a[j-1] );
//
//  ...and if the return value of INDX is 0, we're done.
//
    }
    else
    {
      break;
    }
  }

  i4vec_print ( N, a, "  Sorted array:" );

  return;

# undef N
}
//****************************************************************************80

void test1805 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1805 tests SIMPLEX_VOLUME_ND and TETRAHEDRON_VOLUME_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double tetra[DIM_NUM*(DIM_NUM+1)] = {
     0.000000,  0.942809, -0.333333,
    -0.816496, -0.816496, -0.333333,
     0.816496, -0.816496, -0.333333,
     0.000000,  0.000000,  1.000000 };
  double volume;

  cout << "\n";
  cout << "TEST1805\n";
  cout << "  For an N-dimensional simplex,\n";
  cout << "  SIMPLEX_VOLUME_ND computes the volume.\n";
  cout << "  Here, we check the routine by comparing it\n";
  cout << "  with TETRAHEDRON_VOLUME_3D.\n";

  r8mat_transpose_print ( DIM_NUM, DIM_NUM+1, tetra, "  Simplex vertices:" );

  volume = tetrahedron_volume_3d ( tetra );

  cout << "\n";
  cout << "  Volume computed by TETRAHEDRON_VOLUME_3D:\n";
  cout << "  " << volume << "\n";

  volume = simplex_volume_nd ( DIM_NUM, tetra );

  cout << "\n";
  cout << "  Volume computed by SIMPLEX_VOLUME_ND:\n";
  cout << "  " << volume << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test181 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST181 tests SPHERE_DIA2IMP_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double p1[DIM_NUM] = { -1.0, -1.0, 4.0 };
  double p2[DIM_NUM] = {  5.0,  7.0, 4.0 };
  double pc[DIM_NUM];
  double r;

  cout << "\n";
  cout << "TEST181\n";
  cout << "  SPHERE_DIA2IMP_3D converts a sphere from\n";
  cout << "    diameter to implicit form.\n";

  r8vec_print ( DIM_NUM, p1, "  Point P1:" );
  r8vec_print ( DIM_NUM, p2, "  Point P2:" );

  sphere_dia2imp_3d ( p1, p2, &r, pc );

  cout << "\n";
  cout << "    Radius: " << r << "\n";

  r8vec_print ( DIM_NUM, pc, "  The center:" );

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test182 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST182 tests SPHERE_EXP_CONTAINS_POINT_3D, SPHERE_IMP_CONTAINS_POINT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 4

  int i;
  bool inside;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
    1.0, 2.0, 3.0,
    7.0, 2.0, 3.0,
    1.0, 5.0, 3.0,
    2.5, 3.5, 4.5 };
  double pc[DIM_NUM] = {  1.0, 2.0, 3.0 };
  double p1[DIM_NUM] = {  4.0, 2.0, 3.0 };
  double p2[DIM_NUM] = {  1.0, 5.0, 3.0 };
  double p3[DIM_NUM] = {  1.0, 2.0, 6.0 };
  double p4[DIM_NUM] = { -2.0, 2.0, 3.0 };
  double r = 3.0;
  int test;

  cout << "\n";
  cout << "TEST182\n";
  cout << "  SPHERE_EXP_CONTAINS_POINT_3D determines if a\n";
  cout << "    point is within an explicit sphere;\n";
  cout << "  SPHERE_IMP_CONTAINS_POINT_3D determines if a\n";
  cout << "    point is within an implicit sphere;\n";
  cout << "\n";
  cout << "  SPHERE_EXP_CONTAINS_POINT_3D:\n";
  cout << "    Inside      P\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    inside = sphere_exp_contains_point_3d ( p1, p2, p3, p4, p );

    cout << "  " << setw(1) << inside;
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p[i];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  SPHERE_IMP_CONTAINS_POINT_3D:\n";
  cout << "    Inside      P\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    inside = sphere_imp_contains_point_3d ( r, pc, p );

    cout << "  " << setw(1) << inside;
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p[i];
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test183 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST183 tests SPHERE_EXP_POINT_NEAR_3D and SPHERE_IMP_POINT_NEAR_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 4

  int i;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
    1.0, 2.0, 3.0,
    7.0, 2.0, 3.0,
    1.0, 5.0, 3.0,
    2.5, 3.5, 4.5 };
  double p1[DIM_NUM] = {  4.0, 2.0, 3.0 };
  double p2[DIM_NUM] = {  1.0, 5.0, 3.0 };
  double p3[DIM_NUM] = {  1.0, 2.0, 6.0 };
  double p4[DIM_NUM] = { -2.0, 2.0, 3.0 };
  double pc[DIM_NUM] = {  1.0, 2.0, 3.0 };
  double pn[DIM_NUM];
  double r = 3.0;
  int test;

  cout << "\n";
  cout << "TEST183\n";
  cout << "  SPHERE_EXP_POINT_NEAR_3D determines if a\n";
  cout << "  point is within an explicit sphere;\n";
  cout << "  SPHERE_IMP_POINT_NEAR_3D determines if a\n";
  cout << "  point is within an implicit sphere;\n";
  cout << "\n";
  cout << "  Sphere radius " << r << "\n";

  r8vec_print ( DIM_NUM, pc, "  Sphere center:" );

  cout << "\n";
  cout << "  SPHERE_EXP_POINT_NEAR_3D:\n";
  cout << "    P          PN\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    sphere_exp_point_near_3d ( p1, p2, p3, p4, p, pn );

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p[i];
    }
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << pn[i];
    }
    cout << "\n";
  }
  cout << "\n";
  cout << "  SPHERE_IMP_POINT_NEAR_3D:\n";
  cout << "    P          PN\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    sphere_imp_point_near_3d ( r, pc, p, pn );

    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p[i];
    }
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << pn[i];
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test1835 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1835 tests SPHERE_EXP2IMP_3D and SPHERE_IMP2EXP_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int i;
  double pc[DIM_NUM] = {  1.0, 2.0, 3.0 };
  double p1[DIM_NUM] = {  4.0, 2.0, 3.0 };
  double p2[DIM_NUM] = {  1.0, 5.0, 3.0 };
  double p3[DIM_NUM] = {  1.0, 2.0, 6.0 };
  double p4[DIM_NUM] = { -2.0, 2.0, 3.0 };
  double r = 3.0;

  cout << "\n";
  cout << "TEST1835\n";
  cout << "  SPHERE_EXP2IMP_3D: explicit sphere => implicit form;\n";
  cout << "  SPHERE_IMP2EXP_3D: implicit sphere => explicit form.\n";

  cout << "\n";
  cout << "  Initial form of explicit sphere:\n";
  cout << "\n";
  cout << "  P1:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p1[i];
  }
  cout << "\n";
  cout << "  P2:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p2[i];
  }
  cout << "\n";
  cout << "  P3:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p3[i];
  }
  cout << "\n";
  cout << "  P4:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p4[i];
  }
  cout << "\n";

  sphere_exp2imp_3d ( p1, p2, p3, p4, &r, pc );

  cout << "\n";
  cout << "  Computed form of implicit sphere:\n";
  cout << "\n";
  cout << "  Imputed radius = " << r << "\n";

  r8vec_print ( DIM_NUM, pc, "  Imputed center:" );

  sphere_imp2exp_3d ( r, pc, p1, p2, p3, p4 );

  cout << "\n";
  cout << "  Computed form of explicit sphere:\n";
  cout << "\n";
  cout << "  P1:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p1[i];
  }
  cout << "\n";
  cout << "  P2:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p2[i];
  }
  cout << "\n";
  cout << "  P3:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p3[i];
  }
  cout << "\n";
  cout << "  P4:";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << p4[i];
  }
  cout << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test1836 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1836 tests SPHERE_EXP2IMP_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  int n = N;
  double p[N*(N+1)] = {
    4.0, 2.0, 3.0, 
    1.0, 5.0, 3.0, 
    1.0, 2.0, 6.0, 
   -2.0, 2.0, 3.0 };
  double pc[N];
  double pc_true[N] = { 1.0, 2.0, 3.0 };
  double r;
  double r_true = 3.0;

  cout << "\n";
  cout << "TEST1836\n";
  cout << "  SPHERE_EXP2IMP_ND: explicit sphere => implicit form;\n";

  r8mat_transpose_print ( n, n + 1, p, "  Initial form of explicit sphere:" );

  sphere_exp2imp_nd ( n, p, r, pc );

  cout << "\n";
  cout << "  Computed form of implicit sphere:\n";
  cout << "\n";
  cout << "  Imputed radius = " << r << "\n";
  cout << "  True radius =    " << r_true << "\n";

  r8vec_print ( n, pc, "  Imputed center" );

  r8vec_print ( n, pc_true, "  True center" );

  return;
# undef N
}
//****************************************************************************80

void test187 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST187 tests SPHERE_IMP_GRIDFACES_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TRI_MAX 1000

  int nlat = 3;
  int nlong = 4;
  int triangle_num;
  int tri[DIM_NUM*TRI_MAX];

  cout << "\n";
  cout << "TEST187\n";
  cout << "  SPHERE_IMP_GRIDFACES_3D computes gridfaces\n";
  cout << "  on a sphere in 3D.\n";
  cout << "\n";
  cout << "  Number of intermediate latitudes is " << nlat << "\n";
  cout << "  Number of longitudes is " << nlong << "\n";

  sphere_imp_gridfaces_3d ( TRI_MAX, nlat, nlong, &triangle_num, tri );

  cout << "\n";
  cout << "  The number of triangles is " << triangle_num << "\n";

  i4mat_transpose_print ( DIM_NUM, triangle_num, tri, "  Triangle vertices:" );

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test188 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST188 tests SPHERE_IMP_POINT_PROJECT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 4

  int i;
  double p_test[DIM_NUM*TEST_NUM] = {
    2.0, 0.0,  0.0,
    0.0, 4.0,  0.0,
    2.0, 4.0, 10.0,
    3.0, 5.0,  0.0 };
  double *p1;
  double p2[DIM_NUM];
  double pc[DIM_NUM] = { 2.0, 4.0, 0.0 };
  double r = 2.0;
  int test;

  cout << "\n";
  cout << "TEST188\n";
  cout << "  SPHERE_IMP_POINT_PROJECT_3D projects a 3D point\n";
  cout << "  onto a sphere.\n";
  cout << "\n";
  cout << "        P1       projection P2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p_test + test * DIM_NUM;

    sphere_imp_point_project_3d ( r, pc, p1, p2 );

    cout << "\n";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p1[i];
    }
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p2[i];
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test189 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST189 tests SPHERE_IMP_AREA_ND and SPHERE_IMP_VOLUME_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  int dim_num;
  double r = 1.0;
  double volume;

  cout << "\n";
  cout << "TEST189\n";
  cout << "  For a implicit sphere in N dimensions:\n";
  cout << "  SPHERE_IMP_AREA_ND computes the area;\n";
  cout << "  SPHERE_IMP_VOLUME_ND computes the volume.\n";
  cout << "\n";
  cout << "  We use a radius of R = " << r << "\n";
  cout << "\n";
  cout << "       DIM_NUM        Area            Volume\n";
  cout << "\n";

  for ( dim_num = 2; dim_num <= 10; dim_num++ )
  {
    area = sphere_imp_area_nd ( dim_num, r );
    volume = sphere_imp_volume_nd ( dim_num, r );
    cout << "  " << setw(6)  << dim_num
         << "  " << setw(14) << area
         << "  " << setw(14) << volume << "\n";
  }

  return;
}
//****************************************************************************80

void test1892 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1892 tests SPHERE01_POLYGON_AREA and SPHERE01_POLYGON_AREA_KARNEY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double area1;
  double area2;
  double area_test[15] = { 
     0.001,
     0.01,
     0.1,
     1.0,
    10.0,
    100.0,
    1000.0,
    1.0E+04,
    1.0E+05,
    1.0E+06,
    1.0E+07, 
    1.0E+08,
    1.0E+09,
    1.0E+10,
    1.0E+11 };
  int i;
  int i_test;
  double *lat;
  double lat_test[45] = {
    -80.0, -79.999999567493945456, -79.999999783746965784,
    -70.0, -69.999998632295765827, -69.999999316147849276,
    -60.0, -60.000003745612234573, -59.999999999999717257,
    -50.0, -49.999999999998054558, -49.999988155333397067,
    -40.0, -39.999999999986302375, -39.999962543873523006,
    -30.0, -30.000118446637603681, -29.999999999905752108, 
    -20.0, -20.000374561081985306, -19.999999999405844361,
    -10.0,  -9.999999997121604111,  -9.998815532668760469,
      0.0,   0.002162530270410797,   0.004325060543902236, 
     10.0,   9.993161263027051324,  10.00683830531098623,
     20.0,  19.999994058480001852,  20.037454636874350123, 
     30.0,  29.931544224229581478,  30.06831450158650766,
     40.0,  39.567495485769567709,  39.78272828426239574,
     50.0,  49.301837019419162475,  50.669079015386089932,
     60.0,  55.676479026604850271,  57.644171397676860716  };
  double *lon;
  double lon_test[45] = {
       0.0,  0.0,                  -0.000002157012112311,
       0.0,  0.0,                  -0.000003463148577445,
       0.0,  0.000004325061035168,  0.000008650121090885,
       0.0,  0.000021277700652005,  0.000010638847704917, 
       0.0,  0.000056459655628232,  0.000028229812328744,
       0.0,  0.000078964535025114,  0.000157928881554143,
       0.0,  0.000230132212430763,  0.000460263329704454,
       0.0,  0.001388803276510279,  0.000694399107048577, 
       0.0,  0.003745612305703688, -0.00000000000355722, 
       0.0,  0.012027136076063851,  0.012027642262720488, 
       0.0,  0.046026330173423388,  0.023018642559700887, 
       0.0,  0.136676192075098179,  0.136864622375766145, 
       0.0,  0.0,                   0.487409388137601056, 
       0.0,  1.816527281102587718,  1.868836724548604625, 
       0.0,  0.0,                   7.01051280511881837  };
  int n;
  int n_test[15] = { 
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3 };
  double pi = 3.141592653589793;
  double r;
  int test;
  int test_num = 15;
 
  cout << "\n";
  cout << "TEST1892\n";
  cout << "  For a polygon on the surface of a unit sphere in 3D,\n";
  cout << "  SPHERE01_POLYGON_AREA and\n";
  cout << "  SPHERE01_POLYGON_AREA_KARNEY compute the area.\n";
  cout << "\n";
  cout << "     I     N  AREA           AREA1          ERROR1";
  cout << "         AREA2          ERROR2\n";
  cout << "\n";

  i_test = 0;
  r = 20000000.0 / pi;

  for ( test = 0; test < test_num; test++ )
  {
    n = n_test[test];
    area = area_test[test] / r / r;
    lat = new double[n];
    lon = new double[n];
    for ( i = 0; i < n; i++ )
    {
      lat[i] = lat_test[i_test] * pi / 180.0;
      lon[i] = lon_test[i_test] * pi / 180.0;
      i_test = i_test + 1;
    }
    area1 = sphere01_polygon_area ( n, lat, lon );
    area2 = sphere01_polygon_area_karney ( n, lat, lon );
    cout << "  " << setw(4) << test
         << "  " << setw(4) << n
         << "  " << setw(13) << area
         << "  " << setw(13) << area1
         << "  " << setw(13) << fabs ( area - area1 )
         << "  " << setw(13) << area2
         << "  " << setw(13) << fabs ( area - area2 ) << "\n";

    delete [] lat;
    delete [] lon;
  }

  return;
}
//****************************************************************************80

void test1895 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1895 tests SPHERE_UNIT_AREA_ND and SPHERE_UNIT_AREA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double area2;
  int dim_num;
  int n_data;

  cout << "\n";
  cout << "TEST1895:\n";
  cout << "  SPHERE_UNIT_AREA_ND evaluates the area of the unit\n";
  cout << "  sphere in N dimensions.\n";
  cout << "  SPHERE_UNIT_AREA_VALUES returns some test values.\n";
  cout << "\n";
  cout << "    DIM_NUM    Exact          Computed\n";
  cout << "               Area           Area\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_area_values ( n_data, dim_num, area );

    if ( n_data == 0 )
    {
      break;
    }

    area2 = sphere_unit_area_nd ( dim_num );

    cout << "  "<< setw(6)  << dim_num
         << "  "<< setw(14) << area
         << "  "<< setw(14) << area2 << "\n";
  }

  return;
}
//****************************************************************************80

void test190 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST190 tests SPHERE_UNIT_SAMPLE_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double average[DIM_NUM];
  double dot_average;
  int i;
  int j;
  int sample_num = 1000;
  int seed = 123456789;
  double *v;
  double *x;

  cout << "\n";
  cout << "TEST190\n";
  cout << "  For the unit sphere in 2 dimensions (the circle):\n";
  cout << "  SPHERE_UNIT_SAMPLE_2D samples;\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( j = 1; j <= 5; j++ )
  {
    x = sphere_unit_sample_2d ( seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << x[i];
    }
    cout << "\n";
    delete [] x;
  }
  cout << "\n";
  cout << "  Number of sample points = " << sample_num << "\n";

  r8vec_zero ( DIM_NUM, average );

  for ( j = 1; j <= sample_num; j++ )
  {
    x = sphere_unit_sample_2d ( seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      average[i] = average[i] + x[i];
    }
    delete [] x;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    average[i] = average[i] / ( double ) ( sample_num );
  }
  cout << "\n";
  cout << "  Now average the points, which should get a value\n";
  cout << "  close to zero, and closer as sample_num increases.\n";
  cout << "\n";
  cout << "  Average:        ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(12) << average[i];
  }
  cout << "\n";

  cout << "\n";
  cout << "  Now choose a random direction, sample the same\n";
  cout << "  number of points, and compute the dot product with\n";
  cout << "  the direction.\n";
  cout << "  Take the absolute value of each dot product \n";
  cout << "  and sum and average.\n";
  cout << "\n";
  cout << "  We expect a value near 2 / PI = 0.6366...\n";

  for ( j = 1; j <= 5; j++ )
  {
    v = sphere_unit_sample_2d ( seed );

    dot_average = 0.0;

    for ( i = 1; i <= sample_num; i++ )
    {
      x = sphere_unit_sample_2d ( seed );
      dot_average = dot_average + fabs ( r8vec_dot_product ( DIM_NUM, x, v ) );
      delete [] x;
    }
    dot_average = dot_average / ( double ) ( sample_num );

    cout << "\n";
    cout << "  V:                ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << v[i];
    }
    cout << "\n";

    cout << "  Average |(XdotV)| " << dot_average << "\n";
    delete [] v;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test191 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST191 tests SPHERE_UNIT_SAMPLE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double average[DIM_NUM];
  double dot_average;
  int i;
  int j;
  int k;
  int sample_num = 1000;
  int seed = 123456789;
  double *v;
  double *x;

  cout << "\n";
  cout << "TEST191\n";
  cout << "  For the unit sphere in 3 dimensions:\n";
  cout << "  SPHERE_UNIT_SAMPLE_3D samples;\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( j = 1; j <= 5; j++ )
  {
    x = sphere_unit_sample_3d ( seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << x[i];
    }
    delete [] x;
  }

  r8vec_zero ( DIM_NUM, average );

  for ( j = 1; j <= sample_num; j++ )
  {
    x = sphere_unit_sample_3d ( seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      average[i] = average[i] + x[i];
    }
    delete [] x;
  }

  for ( i = 0; i < DIM_NUM; i++ )
  {
    average[i] = average[i] / ( double ) ( sample_num );
  }
  cout << "\n";
  cout << "  Now average the points, which should get a value\n";
  cout << "  close to zero, and closer as sample_num increases.\n";
  cout << "\n";
  cout << "  Average:        ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << average[i];
  }
  cout << "\n";

  cout << "\n";
  cout << "  Now choose a random direction, sample the same\n";
  cout << "  number of points, and compute the dot product with\n";
  cout << "  the direction.\n";
  cout << "  Take the absolute value of each dot product \n";
  cout << "  and sum and average.\n";

  for ( k = 1; k <= 5; k++ )
  {
    v = sphere_unit_sample_3d ( seed );

    dot_average = 0.0;

    for ( j = 1; j <= sample_num; j++ )
    {
      x = sphere_unit_sample_3d ( seed );
      dot_average = dot_average + fabs ( r8vec_dot_product ( DIM_NUM, x, v ) );
      delete [] x;
    }

    dot_average = dot_average / ( double ) ( sample_num );

    cout << "\n";
    cout << "  V:                ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << v[i];
    }
    cout << "\n";
    cout << "  Average |(XdotV)| " << dot_average << "\n";
    delete [] v;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test192 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST192 tests SPHERE_UNIT_SAMPLE_3D_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double average[DIM_NUM];
  double dot_average;
  int i;
  int j;
  int k;
  int sample_num = 1000;
  int seed = 123456789;
  double *v;
  double *x;

  cout << "\n";
  cout << "TEST192\n";
  cout << "  For the unit sphere in 3 dimensions:\n";
  cout << "  SPHERE_UNIT_SAMPLE_3D_2 samples;\n";
  cout << "\n";
  cout << "  Warning: SPHERE_UNIT_SAMPLE_3D_2 is NOT a good code!\n";
  cout << "  I only implemented it for comparison.\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( j = 1; j <= 5; j++ )
  {
    x = sphere_unit_sample_3d_2 ( seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << x[i];
    }
    delete [] x;
  }

  r8vec_zero ( DIM_NUM, average );

  for ( j = 1; j <= sample_num; j++ )
  {
    x = sphere_unit_sample_3d_2 ( seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      average[i] = average[i] + x[i];
    }
    delete [] x;
  }

  for ( i = 0; i < DIM_NUM; i++ )
  {
    average[i] = average[i] / ( double ) ( sample_num );
  }
  cout << "\n";
  cout << "  Now average the points, which should get a value\n";
  cout << "  close to zero, and closer as sample_num increases.\n";
  cout << "\n";
  cout << "  Average:        ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << average[i];
  }
  cout << "\n";

  cout << "\n";
  cout << "  Now choose a random direction, sample the same\n";
  cout << "  number of points, and compute the dot product with\n";
  cout << "  the direction.\n";
  cout << "  Take the absolute value of each dot product \n";
  cout << "  and sum and average.\n";

  for ( k = 1; k <= 5; k++ )
  {
    v = sphere_unit_sample_3d_2 ( seed );

    dot_average = 0.0;

    for ( j = 1; j <= sample_num; j++ )
    {
      x = sphere_unit_sample_3d_2 ( seed );
      dot_average = dot_average + fabs ( r8vec_dot_product ( DIM_NUM, x, v ) );
      delete [] x;
    }

    dot_average = dot_average / ( double ) ( sample_num );

    cout << "\n";
    cout << "  V:                ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << v[i];
    }
    cout << "\n";
    cout << "  Average |(XdotV)| " << dot_average << "\n";
    delete [] v;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test193 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST193 tests SPHERE_UNIT_SAMPLE_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double average[DIM_NUM];
  double dot_average;
  int i;
  int j;
  int k;
  int sample_num = 1000;
  int seed = 123456789;
  double *v;
  double *x;

  cout << "\n";
  cout << "TEST193\n";
  cout << "  For the unit sphere in N dimensions:\n";
  cout << "  SPHERE_UNIT_SAMPLE_ND samples;\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( j = 1; j <= 5; j++ )
  {
    x = sphere_unit_sample_nd ( DIM_NUM, seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << x[i];
    }
    cout << "\n";
    delete [] x;
  }

  cout << "\n";
  cout << "  Spatial dimension = " << DIM_NUM << "\n";
  cout << "  Number of sample points = " << sample_num << "\n";

  r8vec_zero ( DIM_NUM, average );

  for ( j = 1; j <= sample_num; j++ )
  {
    x = sphere_unit_sample_nd ( DIM_NUM, seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      average[i] = average[i] + x[i];
    }
    delete [] x;
  }

  for ( i = 0; i < DIM_NUM; i++ )
  {
    average[i] = average[i] / ( double ) ( sample_num );
  }
  cout << "\n";
  cout << "  Now average the points, which should get a value\n";
  cout << "  close to zero, and closer as N increases.\n";
  cout << "\n";
  cout << "  Average:        \n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << average[i];
  }
  cout << "\n";
  cout << "\n";
  cout << "  Now choose a random direction, sample the same\n";
  cout << "  number of points, and compute the dot product with\n";
  cout << "  the direction.\n";
  cout << "  Take the absolute value of each dot product \n";
  cout << "  and sum and average.\n";

  for ( k = 1; k <= 5; k++ )
  {
    v = sphere_unit_sample_nd ( DIM_NUM, seed );

    dot_average = 0.0;

    for ( j = 1; j <= sample_num; j++ )
    {
      x = sphere_unit_sample_nd ( DIM_NUM, seed );
      dot_average = dot_average + fabs ( r8vec_dot_product ( DIM_NUM, x, v ) );
      delete [] x;
    }

    dot_average = dot_average / ( double ) ( sample_num );

    cout << "\n";
    cout << "  V:                ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << v[i];
    }
    cout << "\n";
    cout << "  Average |(XdotV)| " << dot_average << "\n";

    delete [] v;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test194 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST194 tests SPHERE_UNIT_SAMPLE_ND_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double average[DIM_NUM];
  double dot_average;
  int i;
  int j;
  int k;
  int sample_num = 1000;
  int seed = 123456789;
  double *v;
  double *x;

  cout << "\n";
  cout << "TEST194\n";
  cout << "  For the unit sphere in N dimensions:\n";
  cout << "  SPHERE_UNIT_SAMPLE_ND_2 samples;\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( j = 1; j <= 5; j++ )
  {
    x = sphere_unit_sample_nd_2 ( DIM_NUM, seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << x[i];
    }
    cout << "\n";
    delete [] x;
  }

  cout << "\n";
  cout << "  Spatial dimension = " << DIM_NUM << "\n";
  cout << "  Number of sample points = " << sample_num << "\n";

  r8vec_zero ( DIM_NUM, average );

  for ( j = 1; j <= sample_num; j++ )
  {
    x = sphere_unit_sample_nd_2 ( DIM_NUM, seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      average[i] = average[i] + x[i];
    }
    delete [] x;
  }

  for ( i = 0; i < DIM_NUM; i++ )
  {
    average[i] = average[i] / ( double ) ( sample_num );
  }
  cout << "\n";
  cout << "  Now average the points, which should get a value\n";
  cout << "  close to zero, and closer as N increases.\n";
  cout << "\n";
  cout << "  Average:        \n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << average[i];
  }
  cout << "\n";
  cout << "\n";
  cout << "  Now choose a random direction, sample the same\n";
  cout << "  number of points, and compute the dot product with\n";
  cout << "  the direction.\n";
  cout << "  Take the absolute value of each dot product \n";
  cout << "  and sum and average.\n";

  for ( k = 1; k <= 5; k++ )
  {
    v = sphere_unit_sample_nd_2 ( DIM_NUM, seed );

    dot_average = 0.0;

    for ( j = 1; j <= sample_num; j++ )
    {
      x = sphere_unit_sample_nd_2 ( DIM_NUM, seed );
      dot_average = dot_average + fabs ( r8vec_dot_product ( DIM_NUM, x, v ) );
      delete [] x;
    }

    dot_average = dot_average / ( double ) ( sample_num );

    cout << "\n";
    cout << "  V:                ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << v[i];
    }
    cout << "\n";
    cout << "  Average |(XdotV)| " << dot_average << "\n";

    delete [] v;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test195 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST195 tests SPHERE_UNIT_SAMPLE_ND_3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double average[DIM_NUM];
  double dot_average;
  int i;
  int j;
  int k;
  int sample_num = 1000;
  int seed = 123456789;
  double *v;
  double *x;

  cout << "\n";
  cout << "TEST195\n";
  cout << "  For the unit sphere in N dimensions:\n";
  cout << "  SPHERE_UNIT_SAMPLE_ND_3 samples;\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( j = 1; j <= 5; j++ )
  {
    x = sphere_unit_sample_nd_3 ( DIM_NUM, seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << x[i];
    }
    cout << "\n";
    delete [] x;
  }

  cout << "\n";
  cout << "  Spatial dimension = " << DIM_NUM << "\n";
  cout << "  Number of sample points = " << sample_num << "\n";

  r8vec_zero ( DIM_NUM, average );

  for ( j = 1; j <= sample_num; j++ )
  {
    x = sphere_unit_sample_nd_3 ( DIM_NUM, seed );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      average[i] = average[i] + x[i];
    }
    delete [] x;
  }

  for ( i = 0; i < DIM_NUM; i++ )
  {
    average[i] = average[i] / ( double ) ( sample_num );
  }
  cout << "\n";
  cout << "  Now average the points, which should get a value\n";
  cout << "  close to zero, and closer as N increases.\n";
  cout << "\n";
  cout << "  Average:        \n";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << "  " << setw(8) << average[i];
  }
  cout << "\n";
  cout << "\n";
  cout << "  Now choose a random direction, sample the same\n";
  cout << "  number of points, and compute the dot product with\n";
  cout << "  the direction.\n";
  cout << "  Take the absolute value of each dot product \n";
  cout << "  and sum and average.\n";

  for ( k = 1; k <= 5; k++ )
  {
    v = sphere_unit_sample_nd_3 ( DIM_NUM, seed );

    dot_average = 0.0;

    for ( j = 1; j <= sample_num; j++ )
    {
      x = sphere_unit_sample_nd_3 ( DIM_NUM, seed );
      dot_average = dot_average + fabs ( r8vec_dot_product ( DIM_NUM, x, v ) );
      delete [] x;
    }

    dot_average = dot_average / ( double ) ( sample_num );

    cout << "\n";
    cout << "  V:                ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(8) << v[i];
    }
    cout << "\n";
    cout << "  Average |(XdotV)| " << dot_average << "\n";

    delete [] v;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test1955 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1955 tests SPHERE_UNIT_VOLUME_ND and SPHERE_UNIT_VOLUME_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  int n_data;
  double volume;
  double volume2;

  cout << "\n";
  cout << "TEST1955:\n";
  cout << "  SPHERE_UNIT_VOLUME_ND evaluates the area of the unit\n";
  cout << "  sphere in N dimensions.\n";
  cout << "  SPHERE_UNIT_VOLUME_VALUES returns some test values.\n";
  cout << "\n";
  cout << "     DIM_NUM    Exact          Computed\n";
  cout << "                Volume         Volume\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_volume_values ( n_data, dim_num, volume );

    if ( n_data == 0 )
    {
      break;
    }

    volume2 = sphere_unit_volume_nd ( dim_num );

    cout << "  " << setw(6)  << dim_num
         << "  " << setw(10) << volume
         << "  " << setw(10) << volume2 << "\n";
  }

  return;
}
//****************************************************************************80

void test196 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST196 tests SHAPE_POINT_DIST_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define SIDE_NUM 4
# define TEST_NUM 9

  double dist;
  int i;
  double *p;
  double p1[DIM_NUM] = { 5.0, 0.0 };
  double pc[DIM_NUM] = { 3.0, 0.0 };
  double p_test[DIM_NUM*TEST_NUM] = {
     3.0,  0.0,
     5.0,  0.0,
     4.0,  0.0,
    10.0,  0.0,
     8.0,  5.0,
     6.0,  6.0,
     1.0,  2.0,
     2.5, -0.5,
     4.0, -1.0 };
  int test;

  cout << "\n";
  cout << "TEST196\n";
  cout << "  For a shape in 2D,\n";
  cout << "  SHAPE_POINT_DIST_2D computes the distance\n";
  cout << "    to a point;\n";
  cout << "\n";
  cout << "  Number of sides:  " << SIDE_NUM << "\n";

  r8vec_print ( DIM_NUM, pc, "  Center of square:" );

  r8vec_print ( DIM_NUM, p1, "  Square vertex #1" );

  cout << "\n";
  cout << "  TEST       X            Y            DIST\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    dist = shape_point_dist_2d ( pc, p1, SIDE_NUM, p );

    cout << "  " << setw(6) << test;
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(12) << p[i];
    }
    cout << "  " << setw(12) << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef SIDE_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test197 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST197 tests SHAPE_POINT_DIST_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define SIDE_NUM 6
# define TEST_NUM 8

  double dist;
  double *p;
  double p1[DIM_NUM] = { 5.0, 0.0 };
  double pc[DIM_NUM] = { 3.0, 0.0 };
  double p_test[DIM_NUM*TEST_NUM] = {
     3.0, 0.0,
     5.0, 0.0,
     4.0, 0.0,
    10.0, 0.0,
     4.0, 1.7320508,
     5.0, 3.4641016,
     3.0, 1.7320508,
     3.0, 0.86602539 };
  int test;

  cout << "\n";
  cout << "TEST197\n";
  cout << "  For a shape in 2D,\n";
  cout << "  SHAPE_POINT_DIST_2D computes the distance\n";
  cout << "    to a point;\n";
  cout << "\n";
  cout << "  Number of sides: " << SIDE_NUM << "\n";

  r8vec_print ( DIM_NUM, pc, "  Center of hexagon:" );

  r8vec_print ( DIM_NUM, p1, "  Hexagon vertex #1" );

  cout << "\n";
  cout << "  TEST         X    Y     DIST\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    dist = shape_point_dist_2d ( pc, p1, SIDE_NUM, p ) ;

    cout << "  " << setw(6) << test
         << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef SIDE_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test198 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST198 tests SHAPE_POINT_NEAR_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define SIDE_NUM 6
# define TEST_NUM 8

  double dist;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
     3.0, 0.0,
     5.0, 0.0,
     4.0, 0.0,
    10.0, 0.0,
     4.0, 1.7320508,
     5.0, 3.4641016,
     3.0, 1.7320508,
     3.0, 0.86602539 };
  double p1[DIM_NUM] = { 5.0, 0.0 };
  double pc[DIM_NUM] = { 3.0, 0.0 };
  double pn[DIM_NUM];
  int test;

  cout << "\n";
  cout << "TEST198\n";
  cout << "  For a shape in 2D,\n";
  cout << "  SHAPE_POINT_NEAR_2D computes the nearest\n";
  cout << "    point to a point;\n";
  cout << "\n";
  cout << "  Number of sides: " << SIDE_NUM << "\n";

  r8vec_print ( DIM_NUM, pc, "  Hexagon center:" );

  r8vec_print ( DIM_NUM, p1, "  Hexagon vertex #1" );

  cout << "\n";
  cout << "  TEST       X            Y              PN     Dist\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    shape_point_near_2d ( pc, p1, SIDE_NUM, p, pn, &dist ) ;

    cout << "  " << setw(6) << test
         << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << pn[0]
         << "  " << setw(10) << pn[1]
         << "  " << setw(10) << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef SIDE_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test199 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST199 tests SHAPE_RAY_INT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define SIDE_NUM 6
# define TEST_NUM 4

  double p1[DIM_NUM] = { 5.0, 0.0 };
  double *pa;
  double pa_test[DIM_NUM*TEST_NUM] = {
    3.0,  0.0,
    3.0,  0.0,
    3.0, -1.0,
    3.0, -1.0 };
  double *pb;
  double pb_test[DIM_NUM*TEST_NUM] = {
    4.0,  0.0,
    3.0,  1.0,
    3.0,  1.0,
    7.0,  5.0 };
  double pc[DIM_NUM] = { 3.0, 0.0 };
  double pint[DIM_NUM];
  int test;

  cout << "\n";
  cout << "TEST199\n";
  cout << "  For a shape in 2D,\n";
  cout << "  SHAPE_RAY_INT_2D computes the intersection of\n";
  cout << "  a shape and a ray whose origin is within\n";
  cout << "  the shape.\n";
  cout << "\n";
  cout << "  Number of sides = " << SIDE_NUM << "\n";

  r8vec_print ( DIM_NUM, pc, "  Hexagon center:" );

  cout << "\n";
  cout << "  Hexagon vertex #1:\n";
  cout << "\n";
  cout << "  " << setw(10) << p1[0]
       << "  " << setw(10) << p1[1] << "\n";

  cout << "\n";
  cout << "  TEST       XA          YA          XB"
       << "          YB          XI          YI\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    pa = pa_test + test * DIM_NUM;
    pb = pb_test + test * DIM_NUM;

    shape_ray_int_2d ( pc, p1, SIDE_NUM, pa, pb, pint );

    cout << "  " << setw(6) << test
         << "  " << setw(10) << pa[0]
         << "  " << setw(10) << pa[1]
         << "  " << setw(10) << pb[0]
         << "  " << setw(10) << pb[1]
         << "  " << setw(10) << pint[0]
         << "  " << setw(10) << pint[1] << "\n";
  }

  return;
# undef DIM_NUM
# undef SIDE_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test200 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST200 tests SPHERE_TRIANGLE_SIDES_TO_ANGLES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2005
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double as;
  double b;
  double bs;
  double c;
  double cs;
  double r = 10.0;

  cout << "\n";
  cout << "TEST200\n";
  cout << "  SPHERE_TRIANGLE_SIDES_TO_ANGLES takes the sides of a\n";
  cout << "  spherical triangle and determines the angles.\n";

  as = 121.0 + ( 15.4 / 60.0 );
  bs = 104.0 + ( 54.7 / 60.0 );
  cs =  65.0 + ( 42.5 / 60.0 );

  as = degrees_to_radians ( as );
  bs = degrees_to_radians ( bs );
  cs = degrees_to_radians ( cs );

  as = r * as;
  bs = r * bs;
  cs = r * cs;
//
//  Get the spherical angles.
//
  sphere_triangle_sides_to_angles ( r, as, bs, cs, a, b, c );

  cout << "\n";
  cout << "  A       = " << a << " (radians)\n";
  a = radians_to_degrees ( a );
  cout << "          = " << a << " ( degrees )\n";
  a = 117.0 + ( 58.0 / 60.0 );
  cout << "  Correct = " << a << " (degrees)\n";

  cout << "\n";
  cout << "  B       = " << b << " (radians)\n";
  b = radians_to_degrees ( b );
  cout << "          = " << b << " ( degrees )\n";
  b = 93.0 + ( 13.8 / 60.0 );
  cout << "  Correct = " << b << " (degrees)\n";

  cout << "\n";
  cout << "  C       = " << c << " (radians)\n";
  c = radians_to_degrees ( c );
  cout << "          = " << c << " ( degrees )\n";
  c = 70.0 + ( 20.6 / 60.0 );
  cout << "  Correct = " << c << " (degrees)\n";

  return;
}
//****************************************************************************80

void test201 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST201 tests STRING_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define VEC_NUM 15

  int i;
  int jstrng;
  int order[VEC_NUM];
  double p1[DIM_NUM*VEC_NUM] = {
    0.0, 0.0,
    3.0, 4.0,
    2.0, 2.0,
    3.0, 2.0,
    2.0, 1.0,
    1.0, 1.0,
    0.0, 5.0,
    1.0, 2.0,
    3.0, 2.0,
    0.0, 0.0,
    5.0, 5.0,
    3.0, 3.0,
    2.0, 4.0,
    7.0, 4.0,
    1.0, 0.0 };
  double p2[DIM_NUM*VEC_NUM] = {
    1.0, 1.0,
    2.0, 4.0,
    1.0, 3.0,
    2.0, 3.0,
    2.0, 2.0,
    1.0, 2.0,
    1.0, 6.0,
    1.0, 3.0,
    3.0, 3.0,
    1.0, 0.0,
    6.0, 6.0,
    3.0, 4.0,
    2.0, 3.0,
    5.0, 5.0,
    2.0, 1.0 };
  int string[VEC_NUM];
  int string_num;

  cout << "\n";
  cout << "TEST201\n";
  cout << "  STRING_2D takes a set of line segments, and\n";
  cout << "    strings them together.\n";
  cout << "\n";
  cout << "  I     P1     P2\n";
  cout << "\n";
  for ( i = 0; i < VEC_NUM; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(10) << p1[0+i*2]
         << "  " << setw(10) << p1[1+i*2]
         << "  " << setw(10) << p2[0+i*2]
         << "  " << setw(10) << p1[1+i*2] << "\n";
  }

  string_2d ( VEC_NUM, p1, p2, &string_num, order, string );

  cout << "\n";
  cout << "  Found " << string_num << " groups of segments.\n";
  cout << "\n";
  cout << "  STRING, ORDER, P1, P2\n";
  cout << "\n";

  jstrng = 1;

  for ( i = 0; i < VEC_NUM; i++ )
  {
    if ( jstrng < string[i] )
    {
      cout << "\n";
      jstrng = jstrng + 1;
    }
    cout << "  " << setw(3)  << string[i]
         << "  " << setw(3)  << order[i]
         << "  " << setw(10) << p1[0+i*2]
         << "  " << setw(10) << p1[1+i*2]
         << "  " << setw(10) << p2[0+i*2]
         << "  " << setw(10) << p2[1+i*2] << "\n";
  }

  return;
# undef DIM_NUM
# undef VEC_NUM
}
//****************************************************************************80

void test202 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST202 tests SUPER_ELLIPSE_POINTS_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 24

  double pc[DIM_NUM] = { 5.0, -2.0 };
  double expo;
  double p[DIM_NUM*N];
  double pi = 3.141592653589793;
  double psi;
  double r1;
  double r2;

  r1 = 3.0;
  r2 = 1.0;
  expo = 1.5;
  psi = pi / 6.0;

  cout << "\n";
  cout << "TEST202\n";
  cout << "  SUPER_ELLIPSE_POINTS_2D returns points on a super ellipse;\n";

  r8vec_print ( DIM_NUM, pc, "  Superellipse center:" );

  cout << "\n";
  cout << "  radii R1 = " << r1 << " R2 = " << r2 << "\n";
  cout << "  exponent EXPO = " << expo << "\n";
  cout << "  and angle PSI = " << psi << "\n";

  super_ellipse_points_2d ( pc, r1, r2, expo, psi, N, p );

  r8mat_transpose_print ( DIM_NUM, N, p, "  Sample points:" );

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test203 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST203 tests TETRAHEDRON_CENTROID_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double *centroid;
  double tetra[DIM_NUM*4] = {
     0.000000,  0.942809, -0.333333,
    -0.816496, -0.816496, -0.333333,
     0.816496, -0.816496, -0.333333,
     0.000000,  0.000000,  1.000000 };

  cout << "\n";
  cout << "TEST203\n";
  cout << "  For a tetrahedron in 3D,\n";
  cout << "  TETRAHEDRON_CENTROID_3D computes the centroid;\n";

  r8mat_transpose_print ( DIM_NUM, 4, tetra, "  Tetrahedron vertices:" );

  centroid = tetrahedron_centroid_3d ( tetra );

  r8vec_print ( DIM_NUM, centroid, "  Centroid:" );

  delete [] centroid;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test2031 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2031 tests TETRAHEDRON_CONTAINS_POINT_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 3

  double *c;
  double c_test[4*TEST_NUM] = {
     0.0, 0.1,  0.2, 0.7,
    -1.3, 2.0,  0.2, 0.1,
     0.8, 0.6, -0.5, 0.1 };
  int i;
  bool inside;
  int j;
  double p[DIM_NUM];
  int test;
  double tetra[DIM_NUM*4] = {
     0.000000,  0.942809, -0.333333,
    -0.816496, -0.816496, -0.333333,
     0.816496, -0.816496, -0.333333,
     0.000000,  0.000000,  1.000000 };

  cout << "\n";
  cout << "TEST2031\n";
  cout << "  For a tetrahedron in 3D,\n";
  cout << "  TETRAHEDRON_CONTAINS_POINT_3D finds if a point \n";
  cout << "  is inside;\n";

  r8mat_transpose_print ( DIM_NUM, 4, tetra, "  Tetrahedron vertices:" );

  cout << "\n";
  cout << "  P     Inside_Tetra?\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    c = c_test + test * 4;

    r8vec_zero ( DIM_NUM, p );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      for ( j = 0; j < 4; j++ )
      {
        p[i] = p[i] + tetra[i+j*DIM_NUM] * c[j];
      }
    }

    inside = tetrahedron_contains_point_3d ( tetra, p );

    cout << "  " << setw(12) << p[0]
         << "  " << setw(12) << p[1]
         << "  " << setw(12) << p[2]
         << "  " << setw(1)  << inside << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2032 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2032 tests TETRAHEDRON_CIRCUMSPHERE_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double pc[DIM_NUM];
  double r;
  double tetra[DIM_NUM*4] = {
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.816496580927726 };

  cout << "\n";
  cout << "TEST2032\n";
  cout << "  For a tetrahedron in 3D,\n";
  cout << "  TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere;\n";

  r8mat_transpose_print ( DIM_NUM, 4, tetra, "  Tetrahedron vertices:" );

  tetrahedron_circumsphere_3d ( tetra, &r, pc );

  r8vec_print ( DIM_NUM, pc, "  Circumsphere center:" );

  cout << "\n";
  cout << "  Circumsphere radius is " << r << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test20321 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20321 tests TETRAHEDRON_EDGE_LENGTH_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double *edge_length;
  double tetra[DIM_NUM*4] = {
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.816496580927726 };

  cout << "\n";
  cout << "TEST20321\n";
  cout << "  For a tetrahedron in 3D,\n";
  cout << "  TETRAHEDRON_EDGE_LENGTH_3D computes the edge lengths;\n";

  r8mat_transpose_print ( DIM_NUM, 4, tetra, "  Tetrahedron vertices:" );

  edge_length = tetrahedron_edge_length_3d ( tetra );

  r8vec_print ( 6, edge_length, "  Edge lengths:" );

  delete [] edge_length;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test20322 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20322 tests TETRAHEDRON_INSPHERE_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double pc[DIM_NUM];
  double r;
  double tetra[DIM_NUM*4] = {
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.816496580927726 };

  cout << "\n";
  cout << "TEST20322\n";
  cout << "  For a tetrahedron in 3D,\n";
  cout << "  TETRAHEDRON_INSPHERE_3D computes the insphere;\n";

  r8mat_transpose_print ( DIM_NUM, 4, tetra, "  Tetrahedron vertices:" );

  tetrahedron_insphere_3d ( tetra, &r, pc );

  r8vec_print ( DIM_NUM, pc, "  Insphere center:" );

  cout << "\n";
  cout << "  Insphere radius is " << r << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test203224 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST203224 tests TETRAHEDRON_LATTICE_LAYER_POINT_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 July 2009
//
//  Author:
//
//    John Burkardt
//
{
  int c[4];
  int i;
  int j;
  int layer;
  bool more;
  int n = 3;
  int v[3];

  cout << "\n";
  cout << "TEST203224\n";
  cout << "  TETRAHEDRON_LATTICE_LAYER_POINT_NEXT returns the next\n";
  cout << "  point in a tetrahedron lattice layer defined by:\n";
  cout << "\n";
  cout << "    C[3] - 1 < X[0]/C[0] + X[1]/C[1] +X[2]/C[2] <= C[3].\n";

  c[0] = 2;
  c[1] = 3;
  c[2] = 4;
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  C =       ";
  for ( i = 0; i < n; i++)
  {
    cout << "  " << setw(4) << c[i];
  }
  cout << "\n";
 
  for ( layer = 0; layer <= 2; layer++ )
  {
    cout << "\n";
    cout << "  Layer " << layer << "\n";
    cout << "\n";

    c[3] = layer;
    more = false;
    i = 0;

    for ( ; ; )
    {
      tetrahedron_lattice_layer_point_next ( c, v, &more );
      if ( !more )
      {
        cout << "  No more.\n";
        break;
      }
      i = i + 1;
      cout << "  " << setw(4) << i;
      for ( j = 0; j < n; j++ )
      {
        cout << "  " << setw(4) << v[j];
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test203225 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST203225 tests TETRAHEDRON_LATTICE_POINT_NEXT.
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
{
# define N 3

  int c[N+1];
  int i;
  int j;
  bool more;
  int n = N;
  int v[N];

  cout << "\n";
  cout << "TEST203225\n";
  cout << "  TETRAHEDRON_LATTICE_POINT_NEXT returns the next lattice\n";
  cout << "  point in a tetrahedron defined by:\n";
  cout << "\n";
  cout << "    0 <= X(1)/C(1) + X(2)/C(2) + X(3)/C(3) <= C(4).\n";

  for ( i = 0; i < n + 1; i++ )
  {
    c[i] = n + 1 - i;
  }
  for ( i = 0; i < n; i++ )
  {
    v[i] = 0;
  }
  more = false;

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  C = ";
  for ( i = 0; i < n + 1; i++ )
  {
    cout << "  " << setw(4) << c[i];
  }
  cout << "\n";
  cout << "\n";

  i = 0;

  for ( ; ; )
  {
    tetrahedron_lattice_point_next ( c, v, &more );
    if ( !more )
    {
      cout << "  No more.\n";
      break;
    }
    i = i + 1;
    cout << "  " << setw(4) << i;
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << v[j];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test20323 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20323 tests TETRAHEDRON_QUALITY1_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  double quality;
  int test;
  double *tetra;
  double tetra_test[DIM_NUM*4*TEST_NUM] = {
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.816496580927726,
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.408248290463863 };

  cout << "\n";
  cout << "TEST20323\n";
  cout << "  For a tetrahedron in 3D,\n";
  cout << "  TETRAHEDRON_QUALITY1_3D computes quality measure #1;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    tetra = tetra_test + test * DIM_NUM * 4;

    r8mat_transpose_print ( DIM_NUM, 4, tetra, "  Tetrahedron vertices:" );

    quality = tetrahedron_quality1_3d ( tetra );

    cout << "\n";
    cout << "  Tetrahedron quality is " << quality << "\n";
 }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test203232 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST203232 tests TETRAHEDRON_QUALITY2_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  double quality2;
  int test;
  double *tetra;
  double tetra_test[DIM_NUM*4*TEST_NUM] = {
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.816496580927726,
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.408248290463863 };

  cout << "\n";
  cout << "TEST203232\n";
  cout << "  For a tetrahedron in 3D,\n";
  cout << "  TETRAHEDRON_QUALITY2_3D computes quality measure #2;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    tetra = tetra_test + test * DIM_NUM * 4;

    r8mat_transpose_print ( DIM_NUM, 4, tetra, "  Tetrahedron vertices:" );

    quality2 = tetrahedron_quality2_3d ( tetra );

    cout << "\n";
    cout << "  Tetrahedron quality is " << quality2 << "\n";
 }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test203233 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST203233 tests TETRAHEDRON_QUALITY3_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  double quality3;
  int test;
  double *tetra;
  double tetra_test[DIM_NUM*4*TEST_NUM] = {
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.816496580927726,
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.408248290463863 };

  cout << "\n";
  cout << "TEST203233\n";
  cout << "  For a tetrahedron in 3D,\n";
  cout << "  TETRAHEDRON_QUALITY3_3D computes quality measure #3;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    tetra = tetra_test + test * DIM_NUM * 4;

    r8mat_transpose_print ( DIM_NUM, 4, tetra, "  Tetrahedron vertices:" );

    quality3 = tetrahedron_quality3_3d ( tetra );

    cout << "\n";
    cout << "  Tetrahedron quality is " << quality3 << "\n";
 }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test203234 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST203234 tests TETRAHEDRON_QUALITY4_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 2

  double quality4;
  int test;
  double *tetra;
  double tetra_test[DIM_NUM*4*TEST_NUM] = {
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.816496580927726,
     0.577350269189626,  0.0, 0.0,
    -0.288675134594813,  0.5, 0.0,
    -0.288675134594813, -0.5, 0.0,
     0.0,                0.0, 0.408248290463863 };

  cout << "\n";
  cout << "TEST203234\n";
  cout << "  For a tetrahedron in 3D,\n";
  cout << "  TETRAHEDRON_QUALITY4_3D computes quality measure #4;\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    tetra = tetra_test + test * DIM_NUM * 4;

    r8mat_transpose_print ( DIM_NUM, 4, tetra, "  Tetrahedron vertices:" );

    quality4 = tetrahedron_quality4_3d ( tetra );

    cout << "\n";
    cout << "  Tetrahedron quality is " << quality4 << "\n";
 }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test203235 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST203235 tests TETRAHEDRON_RHOMBIC_SIZE_3D, TETRAHEDRON_RHOMBIC_SHAPE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int edge_num;
  int face_num;
  int *face_order;
  int face_order_max;
  int *face_point;
  int point_num;
  double *point_coord;

  cout << "\n";
  cout << "TEST203235\n";
  cout << "  For the rhombic tetrahedron,\n";
  cout << "  TETRAHEDRON_RHOMBIC_SIZE_3D returns dimension information;\n";
  cout << "  TETRAHEDRON_RHOMBIC_SHAPE_3D returns face and order information.\n";
  cout << "  SHAPE_PRINT_3D prints this information.\n";
//
//  Get the data sizes.
//
  tetrahedron_rhombic_size_3d ( &point_num, &edge_num, &face_num, 
    &face_order_max );

  cout << "\n";
  cout << "    Number of vertices: " << point_num << "\n";
  cout << "    Number of edges   : " << edge_num << "\n";
  cout << "    Number of faces   : " << face_num << "\n";
  cout << "    Maximum face order: " << face_order_max << "\n";
//
//  Make room for the data.
//
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];
  point_coord = new double[DIM_NUM*point_num];
//
//  Get the data.
//
  tetrahedron_rhombic_shape_3d ( point_num, face_num, face_order_max,
    point_coord, face_order, face_point );
//
//  Print the data.
//
  shape_print_3d ( point_num, face_num, face_order_max,
    point_coord, face_order, face_point );

  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test20324 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20324 tests TETRAHEDRON_SAMPLE_3D, TETRAHEDRON_BARYCENTRIC_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 10

  double p[DIM_NUM];
  int seed = 123456789;
  double t[DIM_NUM*4] = {
     1.0, 4.0, 3.0,
     2.0, 4.0, 3.0,
     1.0, 6.0, 3.0,
     1.0, 4.0, 4.0 };
  int test;
  double *xsi;

  cout << "\n";
  cout << "TEST20324\n";
  cout << "  TETRAHEDRON_SAMPLE_3D samples a tetrahedron.\n";
  cout << "  TETRAHEDRON_BARYCENTRIC_3D converts Cartesian to\n";
  cout << "  barycentric coordinates.\n";
  cout << "\n";
  cout << "  We are computing the barycentric coordinates just to\n";
  cout << "  verify that the points are inside the tetrahedron.\n";

  r8mat_transpose_print ( DIM_NUM, DIM_NUM+1, t, "  Tetrahedron vertices" );

  cout << "\n";
  cout << "     P      Barycentric:\n";
  cout << "\n";

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    tetrahedron_sample_3d ( t, 1, seed, p );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  " << setw(8) << p[2];

    xsi = tetrahedron_barycentric_3d ( t, p );

    cout << "    "
         << "  " << setw(8) << xsi[0]
         << "  " << setw(8) << xsi[1]
         << "  " << setw(8) << xsi[2]
         << "  " << setw(8) << xsi[3] << "\n";

    delete [] xsi;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test20325 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20325 tests TETRAHEDRON_SIZE_3D, TETRAHEDRON_SHAPE_3D, SHAPE_PRINT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int edge_num;
  int face_num;
  int *face_order;
  int face_order_max;
  int *face_point;
  int point_num;
  double *point_coord;

  cout << "\n";
  cout << "TEST20325\n";
  cout << "  For the tetrahedron,\n";
  cout << "  TETRAHEDRON_SIZE_3D returns dimension information;\n";
  cout << "  TETRAHEDRON_SHAPE_3D returns face and order information.\n";
  cout << "  SHAPE_PRINT_3D prints this information.\n";
//
//  Get the data sizes.
//
  tetrahedron_size_3d ( &point_num, &edge_num, &face_num, &face_order_max );

  cout << "\n";
  cout << "    Number of vertices: " << point_num << "\n";
  cout << "    Number of edges   : " << edge_num << "\n";
  cout << "    Number of faces   : " << face_num << "\n";
  cout << "    Maximum face order: " << face_order_max << "\n";
//
//  Make room for the data.
//
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];
  point_coord = new double[DIM_NUM*point_num];
//
//  Get the data.
//
  tetrahedron_shape_3d ( point_num, face_num, face_order_max, point_coord,
    face_order, face_point );
//
//  Print the data.
//
  shape_print_3d ( point_num, face_num, face_order_max,
    point_coord, face_order, face_point );

  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void tetrahedron_solid_angles_3d_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_SOLID_ANGLES_3D tests TETRAHEDRON_VOLUME_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *angle;
  double t1[3*4] = {
     0.000000,  0.942809, -0.333333,
    -0.816496, -0.816496, -0.333333,
     0.816496, -0.816496, -0.333333,
     0.000000,  0.000000,  1.000000 };
  double t2[3*4] = {
     0.000000,  0.000000,  0.000000, 
     1.000000,  0.000000,  0.000000, 
     0.000000,  1.000000,  0.000000, 
     0.000000,  0.000000,  1.000000 };
  double t3[3*4] = {
     0.000000,  0.000000,  0.000000, 
     1.000000,  0.000000,  0.000000, 
     0.000000,  2.000000,  0.000000, 
     0.000000,  0.000000,  4.000000 };
  double t4[3*4] = {
     0.000000,  0.000000,  0.000000, 
     1.000000,  0.000000,  0.000000, 
     0.000000,  1.000000,  0.000000, 
     1.000000,  1.000000,  1.000000 };

  cout << "\n";
  cout << "TETRAHEDRON_SOLID_ANGLES_3D_TEST\n";
  cout << "  TETRAHEDRON_SOLID_ANGLES_3D computes the solid angles\n";
  cout << "  associated with the vertices of a tetrahedron in 3D.\n";

  r8mat_transpose_print ( 3, 4, t1, "  Tetrahedron #1" );
  angle = tetrahedron_solid_angles_3d ( t1 );
  r8vec_print ( 4, angle, "  Solid angles for tetrahedron #1" );
  delete [] angle;

  r8mat_transpose_print ( 3, 4, t2, "  Tetrahedron #2" );
  angle = tetrahedron_solid_angles_3d ( t2 );
  r8vec_print ( 4, angle, "  Solid angles for tetrahedron #2" );
  delete [] angle;

  r8mat_transpose_print ( 3, 4, t3, "  Tetrahedron #3" );
  angle = tetrahedron_solid_angles_3d ( t3 );
  r8vec_print ( 4, angle, "  Solid angles for tetrahedron #3" );
  delete [] angle;

  r8mat_transpose_print ( 3, 4, t4, "  Tetrahedron #4" );
  angle = tetrahedron_solid_angles_3d ( t4 );
  r8vec_print ( 4, angle, "  Solid angles for tetrahedron #4" );
  delete [] angle;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test2033 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2033 tests TETRAHEDRON_VOLUME_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double tetra[DIM_NUM*(DIM_NUM+1)] = {
     0.000000,  0.942809, -0.333333,
    -0.816496, -0.816496, -0.333333,
     0.816496, -0.816496, -0.333333,
     0.000000,  0.000000,  1.000000 };
  double volume;

  cout << "\n";
  cout << "TEST2033\n";
  cout << "  For a tetrahedron in 3D,\n";
  cout << "  TETRAHEDRON_VOLUME_3D computes the volume;\n";

  r8mat_transpose_print ( DIM_NUM, DIM_NUM+1, tetra, "  Tetrahedron vertices" );

  volume = tetrahedron_volume_3d ( tetra );

  cout << "\n";
  cout << "  Volume = " << volume << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test204 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST204 tests TMAT_INIT, TMAT_ROT_AXIS, TMAT_ROT_VECTOR, TMAT_SCALE, TMAT_SHEAR, TMAT_TRANS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a[4*4];
  double angle;
  char axis1;
  string axis2;
  double b[4*4];
  int i;
  double s;
  double v[3];

  cout << "\n";
  cout << "TEST204\n";
  cout << "  TMAT geometric transformation matrix routines:\n";
  cout << "  TMAT_INIT initializes,\n";
  cout << "  TMAT_ROT_AXIS for rotation about an axis,\n";
  cout << "  TMAT_ROT_VECTOR for rotation about a vector,\n";
  cout << "  TMAT_SCALE for scaling,\n";
  cout << "  TMAT_SHEAR for shear,\n";
  cout << "  TMAT_TRANS for translation.\n";
//
//  Initialization.
//
  tmat_init ( a );

  cout << "\n";
  cout << "  Initial transformation matrix:\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    cout << "  " << setw(10) << a[i+0*4]
         << "  " << setw(10) << a[i+1*4]
         << "  " << setw(10) << a[i+2*4]
         << "  " << setw(10) << a[i+3*4] << "\n";
  }
//
//  Rotation about an axis.
//
  angle = 30.0;
  axis1 = 'x';
  tmat_rot_axis ( a, b, angle, axis1 );

  cout << "\n";
  cout << "  Transformation matrix for\n";
  cout << "  rotation about " << axis1 << " by " << angle << " degrees.\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    cout << "  " << setw(10) << b[i+0*4]
         << "  " << setw(10) << b[i+1*4]
         << "  " << setw(10) << b[i+2*4]
         << "  " << setw(10) << b[i+3*4] << "\n";
  }

//
//  Rotation about a vector.
//
  angle = 30.0;
  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;
  tmat_rot_vector ( a, b, angle, v );

  cout << "\n";
  cout << "  Transformation matrix for\n";
  cout << "  rotation about " << v[0] << "  " << v[1] << "  " << v[2]
    << " by " << angle << " degrees.\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    cout << "  " << setw(10) << b[i+0*4]
         << "  " << setw(10) << b[i+1*4]
         << "  " << setw(10) << b[i+2*4]
         << "  " << setw(10) << b[i+3*4] << "\n";
  }
//
//  Scaling.
//
  v[0] = 2.0;
  v[1] = 0.5;
  v[2] = 10.0;
  tmat_scale ( a, b, v );

  cout << "\n";
  cout << "  Transformation matrix for\n";
  cout << "  scaling by " << v[0] << "  " << v[1] << "  " << v[2] << "\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    cout << "  " << setw(10) << b[i+0*4]
         << "  " << setw(10) << b[i+1*4]
         << "  " << setw(10) << b[i+2*4]
         << "  " << setw(10) << b[i+3*4] << "\n";
  }
//
//  Shear.
//
  axis2 = "xy";
  s = 0.5;
  tmat_shear ( a, b, axis2, s );

  cout << "\n";
  cout << "  Transformation matrix for\n";
  cout << "  " << axis2 << " shear coefficient of " << s << "\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    cout << "  " << setw(10) << b[i+0*4]
         << "  " << setw(10) << b[i+1*4]
         << "  " << setw(10) << b[i+2*4]
         << "  " << setw(10) << b[i+3*4] << "\n";
  }
//
//  Translation.
//
  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;
  tmat_trans ( a, b, v );

  cout << "\n";
  cout << "  Transformation matrix for\n";
  cout << "  translation by " << v[0] << "  " << v[1] << "  " << v[2] << "\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    cout << "  " << setw(10) << b[i+0*4]
         << "  " << setw(10) << b[i+1*4]
         << "  " << setw(10) << b[i+2*4]
         << "  " << setw(10) << b[i+3*4] << "\n";
  }

  return;
}
//****************************************************************************80

void test205 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST205 tests TMAT_MXP2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define N 4

  double a[4*4];
  double angle;
  char axis1;
  string axis2;
  double b[4*4];
  double point[DIM_NUM*N];
  double point2[DIM_NUM*N];
  double s;
  double v[3];

  cout << "\n";
  cout << "TEST205\n";
  cout << "  TMAT_MXP2 applies a geometric transformation\n";
  cout << "  matrix to a set of points.\n";
//
//  Initialization.
//
  point[0+0*3] = 1.0;
  point[1+0*3] = 0.0;
  point[2+0*3] = 0.0;

  point[0+1*3] = 0.0;
  point[1+1*3] = 1.0;
  point[2+1*3] = 0.0;

  point[0+2*3] = 0.0;
  point[1+2*3] = 0.0;
  point[2+2*3] = 1.0;

  point[0+3*3] = 1.0;
  point[1+3*3] = 1.0;
  point[2+3*3] = 1.0;

  r8mat_transpose_print ( DIM_NUM, N, point, "  Points:" );
//
//  Initialization of transformation matrix.
//
  tmat_init ( a );
//
//  Rotation about an axis.
//
  angle = 30.0;
  axis1 = 'x';
  tmat_rot_axis ( a, b, angle, axis1 );

  tmat_mxp2 ( b, point, point2, N );

  cout << "\n";
  cout << "  Rotation about " << axis1 << " by " << angle << " degrees.\n";

  r8mat_transpose_print ( DIM_NUM, N, point2, "  Transformed points:" );
//
//  Rotation about a vector.
//
  angle = 30.0;
  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;
  tmat_rot_vector ( a, b, angle, v );

  tmat_mxp2 ( b, point, point2, N );

  cout << "\n";
  cout << "  Rotation about " << v[0]
       << "  " << v[1]
       << "  " << v[2]
       << " by " << angle << " degrees.\n";
  r8mat_transpose_print ( DIM_NUM, N, point2, "  Transformed points:" );
//
//  Scaling.
//
  v[0] = 2.0;
  v[1] = 0.5;
  v[2] = 10.0;
  tmat_scale ( a, b, v );

  tmat_mxp2 ( b, point, point2, N );

  cout << "\n";
  cout << "  Scaling by " << v[0]
       << "  " << v[1]
       << "  " << v[2] << "\n";
  r8mat_transpose_print ( DIM_NUM, N, point2, "  Transformed points:" );
//
//  Shear.
//
  axis2 = "xy";
  s = 0.5;
  tmat_shear ( a, b, axis2, s );

  tmat_mxp2 ( b, point, point2, N );

  cout << "\n";
  cout << "  "  << axis2 << " shear coefficient of " << s << ":\n";
  r8mat_transpose_print ( DIM_NUM, N, point2, "  Transformed points:" );
//
//  Translation.
//
  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;
  tmat_trans ( a, b, v );

  tmat_mxp2 ( b, point, point2, N );

  cout << "\n";
  cout << "  Translation by " << v[0]
       << "  " << v[1]
       << "  " << v[2] << "\n";
  r8mat_transpose_print ( DIM_NUM, N, point2, "  Transformed points:" );

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test206 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST206 tests TRIANGLE_ANGLES_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double angle[3];
  int i;
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };

  cout << "\n";
  cout << "TEST206\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_ANGLES_2D computes the angles;\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  triangle_angles_2d ( t, angle );

  cout << "\n";
  cout << "      Radians      Degrees\n";
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << setw(12) << angle[i]
         << "  " << setw(12) << radians_to_degrees ( angle[i] ) << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test20605 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20605 tests TRIANGLE_ANGLES_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double angle[3];
  int i;
  double t[DIM_NUM*3] = {
    1.0,       2.0,       3.0,
    2.4142137, 3.4142137, 3.0,
    1.7071068, 2.7071068, 4.0 };

  cout << "\n";
  cout << "TEST206\n";
  cout << "  For a triangle in 3D,\n";
  cout << "  TRIANGLE_ANGLES_3D computes the angles;\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  triangle_angles_3d ( t, angle );

  cout << "\n";
  cout << "      Radians      Degrees\n";
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << setw(12) << angle[i]
         << "  " << setw(12) << radians_to_degrees ( angle[i] ) << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test2061 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2061 tests TRIANGLE_AREA_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double area;
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };

  cout << "\n";
  cout << "TEST2061\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_AREA_2D computes the area;\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  area = triangle_area_2d ( t );

  cout << "  Area is " << area << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test2062 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2062 tests TRIANGLE_AREA_HERON;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double area;
  int i;
  int j;
  int jp1;
  double s[3];
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };

  cout << "\n";
  cout << "TEST2062\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_AREA_HERON computes the area;\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  for ( j = 0; j < 3; j++ )
  {
    s[j] = 0.0;

    jp1 = ( j + 1 ) % 3;

    for ( i = 0; i < DIM_NUM; i++ )
    {
      s[j] = s[j] + pow ( t[i+j*DIM_NUM] - t[i+jp1*DIM_NUM], 2 );
    }
    s[j] = sqrt ( s[j] );
  }

  r8vec_print ( 3, s, "  Side lengths:" );

  area = triangle_area_heron ( s );

  cout << "  Area is " << area << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test209 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST209 tests TRIANGLE_AREA_3D, TRIANGLE_AREA_3D_2 and TRIANGLE_AREA_3D_3;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double area;
  double t[DIM_NUM*3] = {
    1.0,       2.0,       3.0,
    2.4142137, 3.4142137, 3.0,
    1.7071068, 2.7071068, 4.0 };

  cout << "\n";
  cout << "TEST209\n";
  cout << "  For a triangle in 3D:\n";
  cout << "  TRIANGLE_AREA_3D   computes the area;\n";
  cout << "  TRIANGLE_AREA_3D_2 computes the area;\n";
  cout << "  TRIANGLE_AREA_3D_3 computes the area;\n";

  r8mat_print ( DIM_NUM, 3, t, "  Triangle (vertices are columns)" );

  area = triangle_area_3d ( t );

  cout << "\n";
  cout << "  Area #1 " << area << "\n";

  area = triangle_area_3d_2 ( t );

  cout << "  Area #2 " << area << "\n";

  area = triangle_area_3d_3 ( t );

  cout << "  Area #3 " << area << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test20655 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20655 tests TRIANGLE_BARYCENTRIC_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 7

  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
     0.25,   0.25,
     0.75,   0.25,
     1.00,   1.00,
    11.00,   0.50,
     0.00,   1.00,
     0.50, -10.00,
     0.60,   0.60 };
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };
  int test;
  double *xsi;

  cout << "\n";
  cout << "TEST20655\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_BARYCENTRIC_2D converts XY coordinates\n";
  cout << "    to barycentric XSI coordinates;\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "          P       XSI\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    xsi = triangle_barycentric_2d ( t, p );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  "
         << "  " << setw(8) << xsi[0]
         << "  " << setw(8) << xsi[1]
         << "  " << setw(8) << xsi[2] << "\n";

    delete [] xsi;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2066 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2066 tests TRIANGLE_CENTROID_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  double *centroid;
  double *t;
  double t_test[DIM_NUM*3*TEST_NUM] = {
         0.0,  0.0,
         1.0,  0.0,
         0.0,  1.0,
         0.0,  0.0,
         1.0,  0.0,
         0.5,  0.86602539,
         0.0,  0.0,
         1.0,  0.0,
         0.5, 10.0,
         0.0,  0.0,
         1.0,  0.0,
        10.0,  2.0 };
  int test;

  cout << "\n";
  cout << "TEST2066\n";
  cout << "  For a triangle in 2D:\n";
  cout << "  TRIANGLE_CENTROID_2D computes the centroid.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    t = t_test + test * DIM_NUM * 3;

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    centroid = triangle_centroid_2d ( t );

    r8vec_print ( DIM_NUM, centroid, "  Centroid:" );

    delete [] centroid;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2094 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2094 tests TRIANGLE_CENTROID_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double *centroid;
  double t[DIM_NUM*3] = {
    1.0,       2.0,       3.0,
    2.4142137, 3.4142137, 3.0,
    1.7071068, 2.7071068, 4.0 };

  cout << "\n";
  cout << "TEST2094\n";
  cout << "  For a triangle in 3D:\n";
  cout << "  TRIANGLE_CENTROID_3D computes the centroid.\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  centroid = triangle_centroid_3d ( t );

  r8vec_print ( DIM_NUM, centroid, "  Centroid:" );

  delete [] centroid;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test2101 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2101 tests tests TRIANGLE_CIRCUMCENTER_2D and others.
//
//  Discussion:
//
//    The functions tested include
//    * TRIANGLE_CIRCUMCENTER_2D;
//    * TRIANGLE_CIRCUMCENTER_2D_2;
//    * TRIANGLE_CIRCUMCENTER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 October 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 2
# define TEST_NUM 4

  int m = M;
  double *pc;
  double *t;
  double t_test[M*3*TEST_NUM] = {
         10.0,  5.0,
         11.0,  5.0,
         10.0,  6.0,
         10.0,  5.0,
         11.0,  5.0,
         10.5,  5.86602539,
         10.0,  5.0,
         11.0,  5.0,
         10.5, 15.0,
         10.0,  5.0,
         11.0,  5.0,
        20.0,   7.0 };
  int test;
  int test_num = TEST_NUM;

  cout << "\n";
  cout << "TEST2101\n";
  cout << "  For a triangle in 2D, the circumcenter can be computed by:\n";
  cout << "  TRIANGLE_CIRCUMCENTER_2D;\n";
  cout << "  TRIANGLE_CIRCUMCENTER_2D_2;\n";
  cout << "  TRIANGLE_CIRCUMCENTER (any dimension);\n";

  for ( test = 0; test < test_num; test++ )
  {
    t = t_test + test * m * 3;

    r8mat_transpose_print ( m, 3, t, "  Triangle vertices:" );

    pc = triangle_circumcenter_2d ( t );
    r8vec_print ( m, pc, "  Circumcenter by TRIANGLE_CIRCUMCENTER_2D:" );
    delete [] pc;

    pc = triangle_circumcenter_2d_2 ( t );
    r8vec_print ( m, pc, "  Circumcenter by TRIANGLE_CIRCUMCENTER_2D_2:" );
    delete [] pc;

    pc = triangle_circumcenter ( m, t );
    r8vec_print ( m, pc, "  Circumcenter by TRIANGLE_CIRCUMCENTER:" );
    delete [] pc;
  }

  return;
# undef M
# undef TEST_NUM
}
//****************************************************************************80

void test21011 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21011 tests TRIANGLE_CIRCUMCENTER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 October 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M1 2
# define TEST_NUM 4

  double *a12;
  int i;
  int j;
  int k;
  int m2;
  double *o1;
  double *o2;
  int m1 = M1;
  double pc1[M1];
  double *pc2;
  int seed;
  double t1[M1*3];
  double *t2;
  double t_test[M1*3*TEST_NUM] = {
         10.0,  5.0, 
         11.0,  5.0, 
         10.0,  6.0, 
         10.0,  5.0, 
         11.0,  5.0, 
         10.5,  5.86602539, 
         10.0,  5.0, 
         11.0,  5.0, 
         10.5, 15.0, 
         10.0,  5.0, 
         11.0,  5.0, 
        20.0,   7.0 };
  int test;
  int test_num = TEST_NUM;

  cout << "\n";
  cout << "TEST21011\n";
  cout << "  For a triangle in M dimensions, the circumenter can be computed by:\n";
  cout << "  TRIANGLE_CIRCUMCENTER;\n";
//
//  Vary the dimension.
//
  for ( m2 = 2; m2 <= 5; m2++ )
  {
    seed = 123456789;

    cout << "\n";
    cout << "  M2 = " << m2 << "\n";

    t2 = new double[m2*3];
//
//  Randomly choose a mapping P2 = O2 + A12 * ( P1 - O1 )
//
    a12 = r8mat_uniform_01_new ( m2, m1, seed );
    o1 = r8vec_uniform_01_new ( m1, seed );
    o2 = r8vec_uniform_01_new ( m2, seed );
//
//  Map each M1-dimensional triangle into M2 space.
//
    for ( test = 0; test < test_num; test++ )
    {
      for ( j = 0; j < 3; j++ )
      {
        for ( i = 0; i < m1; i++ )
        {
          t1[i+j*m1] = t_test[i+j*m1+test*m1*3];
        }
      }
      for ( j = 0; j < 3; j++ )
      {
        t1[i+j*m1] = t1[i+j*m1] - o1[i];
      }

      for ( j = 0; j < 3; j++ )
      {
        for ( i = 0; i < m2; i++ )
        {
          t2[i+j*m2] = 0.0;
          for ( k = 0; k < m1; k++ )
          {
            t2[i+j*m2] = t2[i+j*m2] + a12[i+k*m2] * t1[k+j*m1];
          }
        }
      }
      for ( j = 0; j < 3; j++ )
      {
        for ( i = 0; i < m2; i++ )
        {
          t2[i+j*m2] = t2[i+j*m2] + o2[i];
        }
      }

      pc2 = triangle_circumcenter ( m2, t2 );

      r8vec_print ( m2, pc2, "  Circumcenter by TRIANGLE_CIRCUMCENTER:" );
      cout << "\n";
      cout << "  Distances from circumcenter to vertices:\n";
      cout << "\n";
      for ( j = 0; j < 3; j++ )
      {
        cout << "  " << r8vec_norm_affine ( m2, pc2, t2+j*m2 ) << "\n";
      }
      delete [] pc2;
    }
    delete [] a12;
    delete [] o1;
    delete [] o2;
    delete [] t2;
  }
  return;
}
//****************************************************************************80

void test2067 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2067 tests TRIANGLE_CIRCUMCIRCLE_2D and TRIANGLE_CIRCUMCIRCLE_2D_2;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  double pc[DIM_NUM];
  double r;
  double *t;
  double t_test[DIM_NUM*3*TEST_NUM] = {
         0.0,  0.0,
         1.0,  0.0,
         0.0,  1.0,
         0.0,  0.0,
         1.0,  0.0,
         0.5,  0.86602539,
         0.0,  0.0,
         1.0,  0.0,
         0.5, 10.0,
         0.0,  0.0,
         1.0,  0.0,
        10.0,  2.0 };
  int test;

  cout << "\n";
  cout << "TEST2067\n";
  cout << "  For a triangle in 2D:\n";
  cout << "  TRIANGLE_CIRCUMCIRCLE_2D computes the circumcenter.\n";
  cout << "  TRIANGLE_CIRCUMCIRCLE_2D_2 computes the circumcenter.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    t = t_test + test * DIM_NUM * 3;

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    triangle_circumcircle_2d ( t, &r, pc );

    r8vec_print ( DIM_NUM, pc, "  Circumcenter" );
    cout << "  Circumradius: " << r << "\n";

    triangle_circumcircle_2d_2 ( t, &r, pc );

    r8vec_print ( DIM_NUM, pc, "  Circumcenter2" );
    cout << "  Circumradius2: " << r << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test21015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21015 tests TRIANGLE_CIRCUMRADIUS_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  double r;
  double *t;
  double t_test[DIM_NUM*3*TEST_NUM] = {
        0.0, 0.0,
        1.0, 0.0,
        0.0, 1.0,
        0.0, 0.0,
        1.0, 0.0,
        0.5, 0.86602539,
        0.0,  0.0,
        1.0,  0.0,
        0.5, 10.0,
         0.0, 0.0,
         1.0, 0.0,
        10.0, 2.0 };
  int test;

  cout << "\n";
  cout << "TEST21015\n";
  cout << "  For a triangle in 2D:\n";
  cout << "  TRIANGLE_CIRCUMRADIUS_2D computes the circumradius.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    t = t_test + test * DIM_NUM * 3;

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    r = triangle_circumradius_2d ( t );

    cout << "  Circumradius: " << r << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2068 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2068 tests TRIANGLE_CONTAINS_LINE_EXP_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  bool inside;
  double p1[DIM_NUM] = { 3.0, 0.0, -7.0 };
  double p2[DIM_NUM] = { 5.0, 1.0, -2.0 };
  double pint[DIM_NUM];
  double t[DIM_NUM*3] = {
    8.0, 4.0, 2.0,
    9.0, 0.0, 5.0,
    2.0, 1.0, 2.0 };

  cout << "\n";
  cout << "TEST2068\n";
  cout << "  TRIANGLE_CONTAINS_LINE_EXP_3D determines whether \n";
  cout << "    a triangle contains an explicit line in 3D.\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  r8vec_print ( DIM_NUM, p1, "  Line point P1:" );
  r8vec_print ( DIM_NUM, p2, "  Line point P2:" );

  triangle_contains_line_exp_3d ( t, p1, p2, &inside, pint );

  if ( inside )
  {
    cout << "\n";
    cout << "  The triangle contains the line.\n";
    r8vec_print ( DIM_NUM, pint, "  Intersection point:" );
  }
  else
  {
    cout << "\n";
    cout << "  The triangle does not contain the line.\n";
    r8vec_print ( DIM_NUM, pint, "  The intersection point:" );
  }

  cout << "\n";
  cout << "  Expected answer:\n";
  cout << "\n";
  cout << "    The triangle contains the line, and\n";
  cout << "    the intersection point is at:\n";
  cout << "      7, 2, 3.\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test2069 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2069 tests TRIANGLE_CONTAINS_LINE_PAR_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int dim;
  bool inside;
  double norm;
  double p0[DIM_NUM] = {
    3.0, 0.0, -7.0 };
  double pd[DIM_NUM] = {
    2.0, 1.0, 5.0 };
  double pint[DIM_NUM];
  double t[DIM_NUM*3] = {
    8.0, 4.0, 2.0,
    9.0, 0.0, 5.0,
    2.0, 1.0, 2.0 };

  cout << "\n";
  cout << "TEST2069\n";
  cout << "  TRIANGLE_CONTAINS_LINE_PAR_3D determines whether \n";
  cout << "    a triangle \"contains\" a parametric line in 3D.\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  norm = 0.0;
  for ( dim = 0; dim < DIM_NUM; dim++ )
  {
    norm = norm + pow ( pd[dim], 2 );
  }
  norm = sqrt ( norm );

  for ( dim = 0; dim < DIM_NUM; dim++ )
  {
    pd[dim] = pd[dim] / norm;
  }

  r8vec_print ( DIM_NUM, p0, "  Parametric base point P0:" );
  r8vec_print ( DIM_NUM, pd, "  Parametric direction PD:" );

  triangle_contains_line_par_3d ( t, p0, pd, &inside, pint );

  if ( inside )
  {
    cout << "\n";
    cout << "  The triangle contains the line.\n";
    r8vec_print ( DIM_NUM, pint, "  Intersection point:" );
  }
  else
  {
    cout << "\n";
    cout << "  The triangle does not contain the line.\n";
    r8vec_print ( DIM_NUM, pint, "  The intersection point:" );
  }

  cout << "\n";
  cout << "  Expected answer:\n";
  cout << "\n";
  cout << "    The triangle contains the line, and\n";
  cout << "    the intersection point is at:\n";
  cout << "      ( 7, 2, 3 ).\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test207 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST207 tests TRIANGLE_CONTAINS_POINT_2D_1, TRIANGLE_CONTAINS_POINT_2D_2, TRIANGLE_CONTAINS_POINT_2D_3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 7

  bool inside1;
  bool inside2;
  bool inside3;
  int j;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
     0.25,   0.25,
     0.75,   0.25,
     1.00,   1.00,
    11.00,   0.50,
     0.00,   1.00,
     0.50, -10.00,
     0.60,   0.60 };
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };
  double t2[DIM_NUM*3];
  int test;

  cout << "\n";
  cout << "TEST207\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_CONTAINS_POINT_2D_1 reports if a point\n";
  cout << "    is inside a triangle (and doesn't care about\n";
  cout << "    the ordering of the vertices);\n";
  cout << "  TRIANGLE_CONTAINS_POINT_2D_2 reports if a point \n";
  cout << "    is inside a triangle (and DOES care about\n";
  cout << "    the ordering of the vertices);\n";
  cout << "  TRIANGLE_CONTAINS_POINT_2D_3 reports if a point\n";
  cout << "    is inside a triangle (and doesn't care about\n";
  cout << "    the ordering of the vertices);\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "         X         Y  In1  In2  In3\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    inside1 = triangle_contains_point_2d_1 ( t, p );
    inside2 = triangle_contains_point_2d_2 ( t, p );
    inside3 = triangle_contains_point_2d_3 ( t, p );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "    " << setw(1) << inside1
         << "    " << setw(1) << inside2
         << "    " << setw(1) << inside3 << "\n";
  }
//
//  Make a copy of the triangle with vertices in reverse order.
//
  cout << "\n";
  cout << "  Repeat the test, but reverse the triangle vertex\n";
  cout << "  ordering.\n";

  for ( j = 0; j < 3; j++ )
  {
    t2[0+j*2] = t[0+(2-j)*2];
    t2[1+j*2] = t[1+(2-j)*2];
  }

  r8mat_transpose_print ( DIM_NUM, 3, t2, "  Triangle vertices (reversed):" );

  cout << "\n";
  cout << "         X         Y  In1  In2  In3\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    inside1 = triangle_contains_point_2d_1 ( t2, p );
    inside2 = triangle_contains_point_2d_2 ( t2, p );
    inside3 = triangle_contains_point_2d_3 ( t2, p );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "    " << setw(1) << inside1
         << "    " << setw(1) << inside2
         << "    " << setw(1) << inside3 << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2075 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2075 tests TRIANGLE_DIAMETER_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double diameter;
  double *t;
  double t_test[DIM_NUM*3*TEST_NUM] = {
     4.0, 2.0,
     1.0, 5.0,
    -2.0, 2.0,
     4.0, 2.0,
     5.0, 4.0,
     6.0, 6.0,
     4.0, 2.0,
     1.0, 5.0,
     4.0, 2.0 };
  int test;

  cout << "\n";
  cout << "TEST2075\n";
  cout << "  TRIANGLE_DIAMETER_2D computes the diameter of \n";
  cout << "    the SMALLEST circle around the triangle.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    t = t_test + test * DIM_NUM * 3;

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    diameter = triangle_diameter_2d ( t );

    cout << "  Diameter =       " << diameter << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test208 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST208 tests TRIANGLE_GRIDPOINTS_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define GRID_MAX 50

  double g[DIM_NUM*GRID_MAX];
  int grid_num;
  int sub_num;
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };

  sub_num = 3;

  cout << "\n";
  cout << "TEST208\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_GRIDPOINTS_2D produces a set of\n";
  cout << "  gridpoints in or on the triangle.\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  triangle_gridpoints_2d ( t, sub_num, GRID_MAX, &grid_num, g );

  cout << "\n";
  cout << "  Number of grid points is " << grid_num << "\n";

  r8mat_print ( DIM_NUM, grid_num, g, "  Grid points: " );

  return;
# undef DIM_NUM
# undef GRID_MAX
}
//****************************************************************************80

void test2102 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2102 tests TRIANGLE_INCENTER_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  double pc[DIM_NUM];
  double *t;
  double t_test[DIM_NUM*3*TEST_NUM] = {
        0.0, 0.0,
        1.0, 0.0,
        0.0, 1.0,
        0.0, 0.0,
        1.0, 0.0,
        0.5, 0.86602539,
        0.0,  0.0,
        1.0,  0.0,
        0.5, 10.0,
         0.0, 0.0,
         1.0, 0.0,
        10.0, 2.0 };
  int test;

  cout << "\n";
  cout << "TEST2102\n";
  cout << "  For a triangle in 2D:\n";
  cout << "  TRIANGLE_INCENTER_2D computes the incenter.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    t = t_test + test * DIM_NUM * 3;

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    triangle_incenter_2d ( t, pc );

    r8vec_print ( DIM_NUM, pc, "  Incenter" );
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2070 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2070 tests TRIANGLE_INCIRCLE_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double pc[DIM_NUM];
  double r;
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };

  cout << "\n";
  cout << "TEST2070\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_INCIRCLE_2D computes the incircle;\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  triangle_incircle_2d ( t, pc, &r );

  cout << "  Incircle center is (" << pc[0] << "  " << pc[1] << ").\n";
  cout << "  Incircle radius is " << r << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test20701 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20701 tests TRIANGLE_INRADIUS_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double r;
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };

  cout << "\n";
  cout << "TEST20701\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_INRADIUS_2D computes the inradius;\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  r = triangle_inradius_2d ( t );

  cout << "  Incircle radius is " << r << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test2104 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2104 tests TRIANGLE_LATTICE_LAYER_POINT_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 July 2009
//
//  Author:
//
//    John Burkardt
//
{
  int c[3];
  int i;
  int j;
  int layer;
  int n = 2;
  bool more;
  int v[2];

  cout << "\n";
  cout << "TEST2104\n";
  cout << "  TRIANGLE_LATTICE_LAYER_POINT_NEXT returns the next\n";
  cout << "  point in a triangle lattice layer defined by:\n";
  cout << "\n";
  cout << "    C[2] - 1 < X[0]/C[0] + X[1]/C[1] <= C[2].\n";

  c[0] = 2;
  c[1] = 3;
  v[0] = 0;
  v[1] = 0;

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  C =       ";
  for ( i = 0; i < n; i++)
  {
    cout << "  " << setw(4) << c[i];
  }
  cout << "\n";

  for ( layer = 0; layer <= 4; layer++ )
  {
    cout << "\n";
    cout << "  Layer " << layer << "\n";
    cout << "\n";

    c[2] = layer;
    more = false;
    i = 0;

    for ( ; ; )
    {
      triangle_lattice_layer_point_next ( c, v, &more );
      if ( !more )
      {
        cout << "  No more.\n";
        break;
      }
      i = i + 1;
      cout << "  " << setw(4) << i;
      for ( j = 0; j < n; j++ )
      {
        cout << "  " << setw(4) << v[j];
      }
      cout << "\n";
    }
  }
  return;
}

//****************************************************************************80

void test2105 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2105 tests TRIANGLE_LATTICE_POINT_NEXT.
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
{
# define N 2

  int c[N+1];
  int i;
  int j;
  bool more;
  int n = N;
  int v[N];

  cout << "\n";
  cout << "TEST2105\n";
  cout << "  TRIANGLE_LATTICE_POINT_NEXT returns the next lattice\n";
  cout << "  point in a triangle defined by:\n";
  cout << "\n";
  cout << "    0 <= X(1)/C(1) + X(2)/C(2) <= C(3).\n";

  for ( i = 0; i < n + 1; i++ )
  {
    c[i] = n + 1 - i;
  }
  for ( i = 0; i < n; i++ )
  {
    v[i] = 0;
  }
  more = false;

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  C = ";
  for ( i = 0; i < n + 1; i++ )
  {
    cout << "  " << setw(4) << c[i];
  }
  cout << "\n";
  cout << "\n";

  i = 0;

  for ( ; ; )
  {
    triangle_lattice_point_next ( c, v, &more );
    if ( !more )
    {
      cout << "  No more.\n";
      break;
    }
    i = i + 1;
    cout << "  " << setw(4) << i;
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(4) << v[j];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test211 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST211 tests TRIANGLE_ORIENTATION_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  int i;
  double *t;
  double t_test[DIM_NUM*3*TEST_NUM] = {
     4.0,  2.0,
     1.0,  5.0,
    -2.0,  2.0,
     1.0,  5.0,
     4.0,  2.0,
     1.0, -1.0,
     1.0,  5.0,
     2.0,  7.0,
     3.0,  9.0,
     1.0,  5.0,
     4.0,  2.0,
     1.0,  5.0 };
  int test;

  cout << "\n";
  cout << "TEST211\n";
  cout << "  TRIANGLE_ORIENTATION_2D determines orientation\n";
  cout << "    of a triangle.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    t = t_test + test * DIM_NUM * 3;

    i = triangle_orientation_2d ( t );

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    if ( i == 0 )
    {
      cout << "  The points are counterclockwise.\n";
    }
    else if ( i == 1 )
    {
      cout << "  The points are clockwise.\n";
    }
    else if ( i == 2 )
    {
      cout << "  The points are colinear.\n";
    }
    else if ( i == 3 )
    {
      cout << "  The points are not distinct.\n";
    }
    else
    {
      cout << "  The return value makes no sense.\n";
    }
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2103 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2103 tests TRIANGLE_ORTHOCENTER_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2009
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  bool flag;
  double pc[DIM_NUM];
  double *t;
  double t_test[DIM_NUM*3*TEST_NUM] = {
     0.0,  0.0,
     1.0,  0.0,
     0.0,  1.0,
     0.0,  0.0,
     1.0,  0.0,
     0.5,  0.86602539,
     0.0,  0.0,
     1.0,  0.0,
     0.5, 10.0,
     0.0,  0.0,
     1.0,  0.0,
    10.0,  2.0 };
  int test;

  cout << "\n";
  cout << "TEST2103\n";
  cout << "  For a triangle in 2D:\n";
  cout << "  TRIANGLE_ORTHOCENTER_2D computes the orthocenter.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    t = t_test + test * DIM_NUM * 3;

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    triangle_orthocenter_2d ( t, pc, &flag );

    r8vec_print ( DIM_NUM, pc, "  Orthocenter" );
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2071 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2071 tests TRIANGLE_POINT_DIST_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 7

  double dist;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
     0.25,   0.25,
     0.75,   0.25,
     1.00,   1.00,
    11.00,   0.50,
     0.00,   1.00,
     0.50, -10.00,
     0.60,   0.60 };
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };
  int test;

  cout << "\n";
  cout << "TEST2071\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_POINT_DIST_2D computes the distance\n";
  cout << "    to a point;\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "           P            DIST\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    dist = triangle_point_dist_2d ( t, p );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  " << setw(8) << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test20715 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20715 tests TRIANGLE_POINT_DIST_SIGNED_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 7

  double dist_signed;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
     0.25,   0.25,
     0.75,   0.25,
     1.00,   1.00,
    11.00,   0.50,
     0.00,   1.00,
     0.50, -10.00,
     0.60,   0.60 };
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };
  int test;

  cout << "\n";
  cout << "TEST20715\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_POINT_DIST_SIGNED_2D computes signed\n";
  cout << "    distance to a point;\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "           P       DIST_SIGNED\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    dist_signed = triangle_point_dist_signed_2d ( t, p );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  " << setw(8) << dist_signed << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2095 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2095 tests TRIANGLE_POINT_DIST_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 3

  double dist;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
    1.0,       2.0,       3.0,
    1.3535534, 2.3535534, 3.0,
    0.0,       0.0,       0.0 };
  double t[DIM_NUM*3] = {
    1.0,       2.0,       3.0,
    2.4142137, 3.4142137, 3.0,
    1.7071068, 2.7071068, 4.0 };
  int test;

  cout << "\n";
  cout << "TEST2095\n";
  cout << "  For a triangle in 3D:\n";
  cout << "  TRIANGLE_POINT_DIST_3D computes the distance\n";
  cout << "  to a point;\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "                   P                          DIST\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    dist = triangle_point_dist_3d ( t, p );

    cout << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << p[2]
         << "  " << setw(12) << dist << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2072 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2072 tests TRIANGLE_POINT_NEAR_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 7

  double dist;
  double *p;
  double p_test[DIM_NUM*TEST_NUM] = {
     0.25,   0.25,
     0.75,   0.25,
     1.00,   1.00,
    11.00,   0.50,
     0.00,   1.00,
     0.50, -10.00,
     0.60,   0.60 };
  double pn[DIM_NUM];
  double t[DIM_NUM*3] = {
    0.0, 1.0,
    0.0, 0.0,
    1.0, 0.0 };
  int test;

  cout << "\n";
  cout << "TEST2072\n";
  cout << "  For a triangle in 2D,\n";
  cout << "  TRIANGLE_POINT_NEAR_2D computes the nearest\n";
  cout << "    point to a point.\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "           P                PN\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test + test * DIM_NUM;

    triangle_point_near_2d ( t, p, pn, &dist );

    cout << "  " << setw(10) << p[0]
         << "  " << setw(10) << p[1]
         << "  " << setw(10) << pn[0]
         << "  " << setw(10) << pn[1] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2115 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2115 tests TRIANGLE_QUALITY_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  double quality;
  double *t;
  double t_test[DIM_NUM*3*TEST_NUM] = {
     0.0,  0.0,
     1.0,  0.0,
     0.0,  1.0,
     0.0,  0.0,
     1.0,  0.0,
     0.5,  0.86602539,
     0.0,  0.0,
     1.0,  0.0,
     0.5, 10.0,
     0.0,  0.0,
     1.0,  0.0,
    10.0,  2.0 };
  int test;

  cout << "\n";
  cout << "TEST2115\n";
  cout << "  For a triangle in 2D:\n";
  cout << "  TRIANGLE_QUALITY_2D computes the quality.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    t = t_test + test * DIM_NUM * 3;

    r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

    quality = triangle_quality_2d ( t );

    cout << "  Quality = " << quality << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test212 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST212 tests TRIANGLE_SAMPLE, TRIANGLE_XY_TO_XSI_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 10

  double p[DIM_NUM];
  int seed = 123456789;
  double t[DIM_NUM*3] = {
     4.0, 2.0,
     1.0, 5.0,
    -2.0, 2.0 };
  int test;
  double xsi[DIM_NUM+1];

  cout << "\n";
  cout << "TEST212\n";
  cout << "  TRIANGLE_SAMPLE samples a triangle.\n";
  cout << "  TRIANGLE_XY_TO_XSI_2D converts XY to XSI coordinates.\n";
  cout << "\n";
  cout << "  We are computing the XSI coordinates just to verify\n";
  cout << "  that the points are inside the triangle.\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "  Sample points (X,Y) and (XSI1,XSI2,XSI3) coordinates:\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    triangle_sample ( t, 1, seed, p );
    triangle_xy_to_xsi_2d ( t, p, xsi );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  " << setw(8) << xsi[0]
         << "  " << setw(8) << xsi[1]
         << "  " << setw(8) << xsi[2] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test213 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST213 tests TRIANGLE_SAMPLE, TRIANGLE_XY_TO_XSI_2D, TRIANGLE_XSI_TO_XY_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 10

  int i;
  int j;
  double p[DIM_NUM];
  double p2[DIM_NUM];
  int seed = 123456789;
  double t[DIM_NUM*3] = {
     4.0, 2.0,
     1.0, 5.0,
    -2.0, 2.0 };
  int test;
  double xsi[DIM_NUM+1];

  cout << "\n";
  cout << "TEST213\n";
  cout << "  TRIANGLE_SAMPLE samples a triangle.\n";
  cout << "  TRIANGLE_XY_TO_XSI_2D converts XY to XSI coordinates.\n";
  cout << "  TRIANGLE_XSI_TO_XY_2D converts XSI to XY coordinates.\n";
  cout << "\n";
  cout << "  We verify that (X,Y) -> (XSI1,XSI2,XSI3) -> (X,Y)\n";
  cout << "  works properly.\n";

  r8mat_transpose_print ( DIM_NUM, 3, t, "  Triangle vertices:" );

  cout << "\n";
  cout << "  Sample points:\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    if ( test == 1 )
    {
      for ( i = 0; i < DIM_NUM; i++ )
      {
        p[i] = 0.0;
        for ( j = 0; j < 3; j++ )
        {
          p[i] = p[i] + t[i+j*DIM_NUM];
        }
        p[i] = p[i] / 3.0;
      }
    }
    else if ( test == 2 )
    {
      p[0] = 3.0;
      p[1] = 0.0;
    }
    else
    {
      triangle_sample ( t, 1, seed, p );
    }

    triangle_xy_to_xsi_2d ( t, p, xsi );

    triangle_xsi_to_xy_2d ( t, xsi, p2 );

    cout << "\n";
    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  " << setw(8) << xsi[0]
         << "  " << setw(8) << xsi[1]
         << "  " << setw(8) << xsi[2] << "\n";
    cout << "  " << setw(8) << p2[0]
         << "  " << setw(8) << p2[1] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test219 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST219 tests TUBE_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  double dist;
  double dist_test[TEST_NUM] = { 0.5, 0.5, 1.0, 1.0 };
  int j;
  int n;
  int nlo;
  int n_test[TEST_NUM] = { 4, 5, 5, 5 };
  double *p;
  double p_test[DIM_NUM*19] = {
     0.0,  0.0,
     4.0,  3.0,
     4.0,  0.0,
     0.0,  0.0,
     0.0,  0.0,
     2.0,  0.0,
     2.0,  1.0,
     0.0,  1.0,
     0.0,  0.0,
    10.0, 20.0,
    20.0, 20.0,
    10.0, 10.0,
    20.0, 10.0,
    10.0, 20.0,
     0.0,  0.0,
    10.0,  0.0,
    10.0, 10.0,
    10.0,  0.0,
     0.0,  0.0 };
  double *p1;
  double *p2;
  int test;

  cout << "\n";
  cout << "TEST219\n";
  cout << "  TUBE_2D computes corners of a tube of radius\n";
  cout << "    DIST surrounding a sequence of points.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = n_test[test];
    dist = dist_test[test];

    nlo = 0;
    for ( j = 0; j < test; j++ )
    {
      nlo = nlo + n_test[j];
    }

    p = p_test + DIM_NUM * nlo;

    cout << "\n";
    cout << "  Test " << test << "\n";
    cout << "  Number of points N = " << n << "\n";
    cout << "  Tube radius DIST = " << dist << "\n";

    r8mat_transpose_print ( DIM_NUM, n, p, "  Points to surround:" );

    p1 = new double[DIM_NUM*n];
    p2 = new double[DIM_NUM*n];

    tube_2d ( dist, n, p, p1, p2 );

    r8mat_transpose_print ( DIM_NUM, n, p1, "  P1:" );
    r8mat_transpose_print ( DIM_NUM, n, p2, "  P2:" );

    delete [] p1;
    delete [] p2;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test220 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST220 tests VECTOR_DIRECTIONS_ND;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 5

  double angle[DIM_NUM];
  double angle_degrees[DIM_NUM];
  int j;
  int test;
  double *v;
  double v_test[DIM_NUM*TEST_NUM] = {
     1.0,        0.0,
     1.7320508,  1.0,
    -1.7320508,  1.0,
    -1.7320508, -1.0,
     1.7320508, -1.0 };

  cout << "\n";
  cout << "TEST220\n";
  cout << "  VECTOR_DIRECTIONS_ND computes the angles\n";
  cout << "  that a vector makes with the axes.\n";
  cout << "\n";
  cout << "     X        Y       AX       AY       AX       AY\n";
  cout << "                     (__Radians___)  (___Degrees___)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    v = v_test + test * DIM_NUM;

    vector_directions_nd ( DIM_NUM, v, angle );

    for ( j = 0; j < DIM_NUM; j++ )
    {
      angle_degrees[j] = radians_to_degrees ( angle[j] );
    }

    cout << "  " << setw(7) << v[0]
         << "  " << setw(7) << v[1]
         << "  " << setw(7) << angle[0]
         << "  " << setw(7) << angle[1]
         << "  " << setw(7) << angle_degrees[0]
         << "  " << setw(7) << angle_degrees[1] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test221 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST221 tests VECTOR_DIRECTIONS_ND;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 3

  double angle[DIM_NUM];
  double angle_degrees[DIM_NUM];
  int j;
  int test;
  double *v;
  double v_test[DIM_NUM*TEST_NUM] = {
    1.0, 0.0, 0.0,
    1.0, 2.0, 3.0,
    0.0, 0.0, 1.0 };

  cout << "\n";
  cout << "TEST221\n";
  cout << "  VECTOR_DIRECTIONS_ND computes the angles\n";
  cout << "  that a vector makes with the axes.\n";
  cout << "\n";
  cout << "    X       Y       Z      AX      AY      AZ ";
  cout << "     AX      AY      AZ   \n";
  cout << "                         (_____Radians_______)";
  cout << " (_______Degrees_______)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    v = v_test + test * DIM_NUM;

    vector_directions_nd ( DIM_NUM, v, angle );

    for ( j = 0; j < DIM_NUM; j++ )
    {
      angle_degrees[j] = radians_to_degrees ( angle[j] );
    }
    cout << "  " << setw(7) << v[0]
         << "  " << setw(7) << v[1]
         << "  " << setw(7) << v[2]
         << "  " << setw(7) << angle[0]
         << "  " << setw(7) << angle[1]
         << "  " << setw(7) << angle[2]
         << "  " << setw(7) << angle_degrees[0]
         << "  " << setw(7) << angle_degrees[1]
         << "  " << setw(7) << angle_degrees[2] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test222 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST222 tests VECTOR_ROTATE_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 3

  double angle;
  double a_test[TEST_NUM] = { 30.0, -45.0, 270.0 };
  int test;
  double *v;
  double v_test[DIM_NUM*TEST_NUM] = {
    1.0, 0.0,
    0.0, 2.0,
    1.0, 1.0 };
  double w[DIM_NUM];

  cout << "\n";
  cout << "TEST222\n";
  cout << "  VECTOR_ROTATE_2D rotates a vector through\n";
  cout << "  a given angle around the origin.\n";
  cout << "\n";
  cout << "    X1      Y1   Angle      X2      Y2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    v = v_test + test * DIM_NUM;

    angle = degrees_to_radians ( a_test[test] );

    vector_rotate_2d ( v, angle, w );

    cout << "  " << setw(7) << v[0]
         << "  " << setw(7) << v[1]
         << "  " << setw(7) << a_test[test]
         << "  " << setw(7) << w[0]
         << "  " << setw(7) << w[1] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2225 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2225 tests VECTOR_ROTATE_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 5

  double axis[DIM_NUM] = { 1.0, 1.0, 1.0 };
  double angle;
  double a_test[TEST_NUM] = { 30.0, -45.0, 90.0, 270.0, 30.0 };
  int test;
  double *v1;
  double v1_test[DIM_NUM*TEST_NUM] = {
    1.0,  0.0,  0.0,
    0.0,  2.0,  0.0,
    0.0,  0.0,  3.0,
    1.0,  1.0,  1.0,
    1.0,  1.0, -2.0 };
  double v2[DIM_NUM];

  cout << "\n";
  cout << "TEST2225\n";
  cout << "  VECTOR_ROTATE_3D rotates a vector through\n";
  cout << "  a given angle around the origin.\n";
  cout << "\n";
  cout << "  Rotations will be about the following axis:\n";
  cout << "\n";
  cout << "  " << setw(8) << axis[0]
       << "  " << setw(8) << axis[1]
       << "  " << setw(8) << axis[2] << "\n";
  cout << "\n";
  cout << "              V1             Angle             V2\n";
  cout << "    ----------------------  ------  ----------------------";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    v1 = v1_test + test * DIM_NUM;

    angle = degrees_to_radians ( a_test[test] );

    vector_rotate_3d ( v1, axis, angle, v2 );

    cout << "  " << setprecision(4) << setw(7) << v1[0]
         << "  " << setprecision(4) << setw(7) << v1[1]
         << "  " << setprecision(4) << setw(7) << v1[2]
         << "  " << setprecision(4) << setw(7) << angle
         << "  " << setprecision(4) << setw(7) << v2[0]
         << "  " << setprecision(4) << setw(7) << v2[1]
         << "  " << setprecision(4) << setw(7) << v2[2] << "\n";
  }
//
//  Test using an axis that is not of unit length!
//
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 2.0;
  cout << "\n";
  cout << "  Rotations will be about the following axis:\n";
  cout << "\n";
  cout << "  " << setw(8) << axis[0]
       << "  " << setw(8) << axis[1]
       << "  " << setw(8) << axis[2] << "\n";
  cout << "\n";
  cout << "              V1             Angle             V2\n";
  cout << "    ----------------------  ------  ----------------------";
  cout << "\n";

  v1[0] = 1.0;
  v1[1] = 1.0;
  v1[2] = 1.0;

  angle = 90.0;
  angle = degrees_to_radians ( angle );

  vector_rotate_3d ( v1, axis, angle, v2 );

  cout << "  " << setprecision(4) << setw(7) << v1[0]
       << "  " << setprecision(4) << setw(7) << v1[1]
       << "  " << setprecision(4) << setw(7) << v1[2]
       << "  " << setprecision(4) << setw(7) << 90.0
       << "  " << setprecision(4) << setw(7) << v2[0]
       << "  " << setprecision(4) << setw(7) << v2[1]
       << "  " << setprecision(4) << setw(7) << v2[2] << "\n";

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test223 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST223 tests VECTOR_ROTATE_BASE_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define TEST_NUM 4

  double angle;
  double a_test[TEST_NUM] = { 30.0, -45.0, 270.0, 20.0 };
  double *p1;
  double p2[DIM_NUM];
  double pb[DIM_NUM] = { 10.0, 5.0 };
  double p_test[DIM_NUM*TEST_NUM] = {
    11.0, 5.0,
    10.0, 7.0,
    11.0, 6.0,
    10.0, 5.0 };
  int test;

  cout << "\n";
  cout << "TEST223\n";
  cout << "  VECTOR_ROTATE_BASE_2D rotates a vector (X1,Y1)\n";
  cout << "  through an angle around a base point (XB,YB).\n";
  cout << "\n";
  cout << "        P1              PB       Angle          P2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p1 = p_test + test * DIM_NUM;

    angle = degrees_to_radians ( a_test[test] );

    vector_rotate_base_2d ( p1, pb, angle, p2 );

    cout << "  " << setw(7) << p1[0]
         << "  " << setw(7) << p1[1]
         << "  " << setw(7) << pb[0]
         << "  " << setw(7) << pb[1]
         << "  " << setw(7) << a_test[test]
         << "  " << setw(7) << p2[0]
         << "  " << setw(7) << p2[1] << "\n";
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test224 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST224 tests VECTOR_SEPARATION_ND;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define TEST_NUM 5

  int test1;
  int test2;
  double theta;
  double theta_deg;
  double *v1;
  double *v2;
  double v_test[DIM_NUM*TEST_NUM] = {
    1.0,  0.0,  0.0,
    1.0,  2.0,  3.0,
    0.0,  0.0,  1.0,
   -3.0,  2.0, -1.0,
   -2.0, -4.0, -6.0 };

  cout << "\n";
  cout << "TEST224\n";
  cout << "  VECTOR_SEPARATION_3D computes the separation angle\n";
  cout << "  between two vectors.\n";
  cout << "\n";
  cout << "    -----Vector 1-----      -----Vector 2-----  \n";
  cout << "   Radians    Degrees\n";
  cout << "\n";

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    v1 = v_test + test1 * DIM_NUM;

    for ( test2 = test1 + 1; test2 < TEST_NUM; test2++ )
    {
      v2 = v_test + test2 * DIM_NUM;

      theta = vector_separation_nd ( DIM_NUM, v1, v2 );

      theta_deg = radians_to_degrees ( theta );

      cout << "  " << setw(7) << v1[0]
           << "  " << setw(7) << v1[1]
           << "  " << setw(7) << v1[2]
           << "  " << setw(7) << v2[0]
           << "  " << setw(7) << v2[1]
           << "  " << setw(7) << v2[2]
           << "  " << setw(7) << theta
           << "  " << setw(7) << theta_deg << "\n";
    }
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test2245 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2245 tests VOXELS_DIST_L1_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int dist;
  int p1[DIM_NUM] = { 1, 1, 5 };
  int p2[DIM_NUM] = { 9, 4, 4 };

  cout << "\n";
  cout << "TEST2245\n";
  cout << "  VOXELS_DIST_L1_ND prints the voxels on a line in 3D.\n";

  cout << "\n";
  cout << "  P1:\n";
  cout << "  " << setw(6) << p1[0]
       << "  " << setw(6) << p1[1]
       << "  " << setw(6) << p1[2] << "\n";


  cout << "\n";
  cout << "  P2:\n";
  cout << "  " << setw(6) << p2[0]
       << "  " << setw(6) << p2[1]
       << "  " << setw(6) << p2[2] << "\n";

  dist = voxels_dist_l1_nd ( DIM_NUM, p1, p2 );

  cout << "\n";
  cout << "  L1 distance = " << dist << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test225 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST225 tests VOXELS_LINE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int n;
  int *v;
  int p1[DIM_NUM] = { 1, 1, 5 };
  int p2[DIM_NUM] = { 9, 4, 4 };

  cout << "\n";
  cout << "TEST225\n";
  cout << "  VOXELS_LINE_3D computes the voxels on a line in 3D\n";
  cout << "  starting at the first voxel, and heading towards\n";
  cout << "  the second one.\n";

  cout << "\n";
  cout << "  Starting voxel:\n";
  cout << "  " << setw(6) << p1[0]
       << "  " << setw(6) << p1[1]
       << "  " << setw(6) << p1[2] << "\n";

  cout << "\n";
  cout << "  Heading voxel:\n";
  cout << "  " << setw(6) << p2[0]
       << "  " << setw(6) << p2[1]
       << "  " << setw(6) << p2[2] << "\n";

  n = voxels_dist_l1_nd ( DIM_NUM, p1, p2 ) + 1;

  cout << "\n";
  cout << "  Number of voxels we will compute is " << n << "\n";

  v = new int[3*n];

  voxels_line_3d ( p1, p2, n, v );

  i4mat_transpose_print ( 3, n, v, "  The voxels:" );

  delete [] v;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test226 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST226 tests VOXELS_REGION_3D.
//
//  Discussion:
//
//    The test region is 8 by 9 by 1 voxels:
//
//    123456789
//  1 .........
//  2 ...11.1..
//  3 ..11111..
//  4 ...11.1..
//  5 ......1..
//  6 .11..11..
//  7 ..1......
//  8 .......1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define LIST_MAX 100
# define NX 8
# define NY 9
# define NZ 1

  int i;
  int ishow[NX*NY*NZ];
  int j;
  int k;
  int l;
  int list[LIST_MAX];
  int list_num;
  int nelements;
  int region;
  int region_num;

  cout << "\n";
  cout << "TEST226\n";
  cout << "  VOXELS_REGION_3D groups voxels into regions.\n";

  i4vec_zero ( NX*NY*NZ, ishow );

  ishow[1+3*NX+0*NX*NY] = 1;
  ishow[1+4*NX+0*NX*NY] = 1;
  ishow[1+6*NX+0*NX*NY] = 1;

  ishow[2+2*NX+0*NX*NY] = 1;
  ishow[2+3*NX+0*NX*NY] = 1;
  ishow[2+4*NX+0*NX*NY] = 1;
  ishow[2+5*NX+0*NX*NY] = 1;
  ishow[2+6*NX+0*NX*NY] = 1;

  ishow[3+3*NX+0*NX*NY] = 1;
  ishow[3+4*NX+0*NX*NY] = 1;
  ishow[3+6*NX+0*NX*NY] = 1;

  ishow[4+6*NX+0*NX*NY] = 1;

  ishow[5+1*NX+0*NX*NY] = 1;
  ishow[5+2*NX+0*NX*NY] = 1;
  ishow[5+5*NX+0*NX*NY] = 1;
  ishow[5+6*NX+0*NX*NY] = 1;

  ishow[6+2*NX+0*NX*NY] = 1;

  ishow[7+7*NX+0*NX*NY] = 1;

  voxels_region_3d ( LIST_MAX, NX, NY, NZ, ishow, &list_num, list,
    &region_num );

  cout << "\n";
  cout << "  Number of regions found = " << region_num << "\n";
  cout << "\n";
  cout << "  The nonzero ISHOW array elements are:\n";
  cout << "\n";

  for ( i = 0; i < NX; i++ )
  {
    for ( j = 0; j < NY; j++ )
    {
      for ( k = 0; k < NZ; k++ )
      {
        l = ishow[i+j*NX+k*NX*NY];
        if ( l != 0 )
        {
          cout << "  " << setw(6) << i+1
               << "  " << setw(6) << j+1
               << "  " << setw(6) << k+1
               << "  " << setw(6) << l << "\n";
        }
      }
    }
  }

  if ( LIST_MAX < list_num )
  {
    cout << "\n";
    cout << "  The stack-based list of regions is unusable.\n";
  }
  else
  {
    cout << "\n";
    cout << "  The stack-based list of regions is:\n";
    cout << "\n";

    region = region_num;

    while ( 0 < list_num )
    {
      nelements = list[list_num-1];
      list_num = list_num - 1;

      cout << "\n";
      cout << "  Region " << region
           << " includes " << nelements << " voxels:\n";
      cout << "\n";

      for ( l = 1; l <= nelements; l++ )
      {
        k = list[list_num-1];
        list_num = list_num - 1;
        j = list[list_num-1];
        list_num = list_num - 1;
        i = list[list_num-1];
        list_num = list_num - 1;
        cout << "  " << setw(6) << i
             << "  " << setw(6) << j
             << "  " << setw(6) << k << "\n";
      }

      region = region - 1;
    }
  }

  return;
# undef DIM_NUM
# undef LIST_MAX
# undef NX
# undef NY
# undef NZ
}
//****************************************************************************80

void test227 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST227 tests VOXELS_STEP_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int i;
  int inc;
  int jnc;
  int knc;
  int v1[DIM_NUM] = { 1, 1, 5 };
  int v2[DIM_NUM];
  int v3[DIM_NUM];

  cout << "\n";
  cout << "TEST227\n";
  cout << "  VOXELS_STEP_3D steps along a line from\n";
  cout << "    one voxel to another.\n";

  i4vec_copy ( DIM_NUM, v1, v2 );

  inc = 7;
  jnc = 3;
  knc = -1;

  cout << "\n";
  cout << "  " << setw(4) << 0
       << "  " << setw(6) << v2[0]
       << "  " << setw(6) << v2[1]
       << "  " << setw(6) << v2[2] << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    voxels_step_3d ( v1, v2, inc, jnc, knc, v3 );
    cout << "  " << setw(4) << i
         << "  " << setw(6) << v3[0]
         << "  " << setw(6) << v3[1]
         << "  " << setw(6) << v3[2] << "\n";
    i4vec_copy ( DIM_NUM, v3, v2 );
  }

  cout << "\n";
  cout << "  Now, as a check, reverse direction and return.\n";
  cout << "\n";

  i4vec_copy ( DIM_NUM, v2, v1 );

  inc = -inc;
  jnc = -jnc;
  knc = -knc;

  i4vec_copy ( DIM_NUM, v1, v2 );

  cout << "  " << setw(4) << 0
       << "  " << setw(6) << v2[0]
       << "  " << setw(6) << v2[1]
       << "  " << setw(6) << v2[2] << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    voxels_step_3d ( v1, v2, inc, jnc, knc, v3 );
    cout << "  " << setw(4) << i
         << "  " << setw(6) << v3[0]
         << "  " << setw(6) << v3[1]
         << "  " << setw(6) << v3[2] << "\n";
    i4vec_copy ( DIM_NUM, v3, v2 );
  }

  return;
# undef DIM_NUM
}
