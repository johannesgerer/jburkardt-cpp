# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "geompack.hpp"

int main ( );
void test005 ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMPACK_PRB tests routines from the GEOMPACK library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "GEOMPACK_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the GEOMPACK library.\n";

  test005 ( );
  test01 ( );
  test02 ( );
  test03 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "GEOMPACK_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests DIAEDG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int node_num = 4;
  int triangle_num = 2;
  int triangle_order = 3;

  double alpha_area;
  double alpha_ave;
  double alpha_min_swapped;
  double alpha_min_unswapped;
  double node_xy[2*4];
  int seed = 123456789;
  bool swap;
  int test;
  int test_num = 10;
  int triangle_node[3*2];
  int value;

  cout << "\n";
  cout << "TEST005\n";
  cout << "  DIAEDG determines whether two triangles\n";
  cout << "  with a common edge need to \"swap\" diagonals.\n";
  cout << "  If swapping is indicated, then ALPHA_MIN should increase.\n";
  cout << "\n";
  cout << "  Swap   ALPHA_MIN   ALPHA_MIN\n";
  cout << "         Unswapped   Swapped\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
//
//  Generate a random quadrilateral (1,2,3,4).
//
    quad_convex_random ( &seed, node_xy );
//
//  Does it need swapping?
//
    value = diaedg (
      node_xy[0+0*2], node_xy[1+0*2],
      node_xy[0+1*2], node_xy[1+1*2],
      node_xy[0+2*2], node_xy[1+2*2],
      node_xy[0+3*2], node_xy[1+3*2] );

    if ( value == 1 )
    {
      swap = false;
    }
    else
    {
      swap = true;
    }
//
//  Compute ALPHA_MIN unswapped.
//
    triangle_node[0+0*3] = 1;
    triangle_node[1+0*3] = 2;
    triangle_node[2+0*3] = 3;
    triangle_node[0+1*3] = 1;
    triangle_node[1+1*3] = 3;
    triangle_node[2+1*3] = 4;

    alpha_measure ( node_num, node_xy, triangle_order, triangle_num,
      triangle_node, &alpha_min_unswapped, &alpha_ave, &alpha_area );
//
//  Compute ALPHA_MIN swapped.
//
    triangle_node[0+0*3] = 1;
    triangle_node[1+0*3] = 2;
    triangle_node[2+0*3] = 4;
    triangle_node[0+1*3] = 2;
    triangle_node[1+1*3] = 3;
    triangle_node[2+1*3] = 4;

    alpha_measure ( node_num, node_xy, triangle_order, triangle_num,
      triangle_node, &alpha_min_swapped, &alpha_ave, &alpha_area );

    if ( false )
    {
      r8mat_transpose_print ( 2, node_num, node_xy, "  Quadrilateral" );
    }

    cout << "     " << swap
         << "  " << setw(10) << alpha_min_unswapped
         << "  " << setw(10) << alpha_min_swapped << "\n";
  }

  return;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests POINTS_DELAUNAY_NAIVE_2D.
//
//  Diagram:
//
//   !....3&11....
//   !............
//   !............
//   X..9.........
//   !.....5......
//   !...........6
//   !.4.2...10...
//   !.....8...12.
//   V............
//   !..7.........
//   !......1.....
//   !............
//   !............
//   !----V----X--
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2009
//
//  Author:
//
//    John Burkardt
//
# define N 12
# define DIM_NUM 2
{
  int i;
  int ntri;
  int *tri;
  double p[DIM_NUM*N] = {
     7.0, 3.0,
     4.0,  7.0,
     5.0, 13.0,
     2.0,  7.0,
     6.0,  9.0,
    12.0, 10.0,
     3.0,  4.0,
     6.0,  6.0,
     3.0, 10.0,
     8.0,  7.0,
     5.0, 13.0,
    10.0,  6.0 };

  cout << "\n";
  cout << "TEST01\n";
  cout << "  POINTS_DELAUNAY_NAIVE_2D computes the Delaunay\n";
  cout << "  triangulation of a set of points.\n";

  r8mat_transpose_print ( DIM_NUM, N, p, "  The points:" );

  tri = points_delaunay_naive_2d ( N, p, &ntri );

  cout << "\n";
  cout << "  Number of triangles is NTRI = " << ntri << "\n";

  i4mat_transpose_print ( 3, ntri, tri, "  The Delaunay triangles:" );

  delete [] tri;

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R82VEC_PART_QUICK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 12
# define DIM_NUM 2

  double *a;
  int i;
  int j;
  int l;
  int r;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  D2VEC_PART_QUICK_A reorders a D2 vector\n";
  cout << "    as part of a quick sort.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8mat_uniform_01_new ( DIM_NUM, N, &seed );

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      a[i+j*DIM_NUM] = 10.0 * a[i+j*DIM_NUM];
    }
  }

  r8mat_transpose_print ( DIM_NUM, N, a, "  Before rearrangment:" );

  r82vec_part_quick_a ( N, a, &l, &r );

  cout << "\n";
  cout << "  Rearranged array\n";
  cout << "  Left index =  " << l << "\n";
  cout << "  Key index =   " << l+1 << "\n";
  cout << "  Right index = " << r << "\n";

  r8mat_transpose_print ( DIM_NUM, l,     a,         "  Left half:" );
  r8mat_transpose_print ( DIM_NUM, 1,     a+2*l,     "  Key:" );
  r8mat_transpose_print ( DIM_NUM, N-l-1, a+2*(l+1), "  Right half:" );

  delete [] a;

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests R82VEC_SORT_QUICK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 12
# define DIM_NUM 2

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int j;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  R82VEC_SORT_QUICK_A sorts a D2 vector\n";
  cout << "    as part of a quick sort.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8mat_uniform_01_new ( DIM_NUM, N, &seed );

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      a[i+j*DIM_NUM] = 10.0 * a[i+j*DIM_NUM];
    }
  }
//
//  For better testing, give a few elements the same first component.
//
  a[0+2*(3-1)] = a[0+2*(5-1)];
  a[0+2*(4-1)] = a[0+2*(12-1)];
//
//  Make two entries equal.
//
  a[0+2*(7-1)] = a[0+2*(11-1)];
  a[1+2*(7-1)] = a[1+2*(11-1)];

  r8mat_transpose_print ( DIM_NUM, N, a, "  Before sorting:" );

  r82vec_sort_quick_a ( N, a );

  r8mat_transpose_print ( DIM_NUM, N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests R8TRIS2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define NODE_NUM 9

  int error;
  double g_xy[DIM_NUM*NODE_NUM] = {
       0.0, 0.0,
       0.0, 1.0,
       0.2, 0.5,
       0.3, 0.6,
       0.4, 0.5,
       0.6, 0.4,
       0.6, 0.5,
       1.0, 0.0,
       1.0, 1.0 };
  int nod_tri[2*NODE_NUM*3];
  int triangle_neighbor[2*NODE_NUM*3];
  int tri_num;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  R8TRIS2 computes the Delaunay triangulation of a\n";
  cout << "  pointset in 2D.\n";
//
//  Set up the Delaunay triangulation.
//
  error = r8tris2 ( NODE_NUM, g_xy, &tri_num, nod_tri, triangle_neighbor );

  if ( error == 0 )
  {
    cout << "\n";
    cout << "  R8TRIS2 computed the Delaunay triangulation with no\n";
    cout << "  errors detected.\n";
  }
  else
  {
    cout << "\n";
    cout << "  R8TRIS2 detected an error condition of index " << error << "\n";
    return;
  }

  triangulation_print ( NODE_NUM, g_xy, tri_num, nod_tri, triangle_neighbor );

  return;
# undef NODE_NUM
# undef DIM_NUM
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests TRIANGLE_CIRCUMCENTER_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double *center;
  int i;
  int ntest = 4;
  double t[DIM_NUM*3];

  cout << "\n";
  cout << "TEST06\n";
  cout << "  For a triangle in 2D:\n";
  cout << "  TRIANGLE_CIRCUMCENTER_2D computes the circumcenter.\n";

  for ( i = 1; i <= ntest; i++ )
  {
    if ( i == 1 )
    {
      t[0+0*2] = 0.0;
      t[1+0*2] = 0.0;
      t[0+1*2] = 1.0;
      t[1+1*2] = 0.0;
      t[0+2*2] = 0.0;
      t[1+2*2] = 1.0;
    }
    else if ( i == 2 )
    {
      t[0+0*2] = 0.0;
      t[1+0*2] = 0.0;
      t[0+1*2] = 1.0;
      t[1+1*2] = 0.0;
      t[0+2*2] = 0.5;
      t[1+2*2] = sqrt ( 3.0 ) / 2.0;
    }
    else if ( i == 3 )
    {
      t[0+0*2] = 0.0;
      t[1+0*2] = 0.0;
      t[0+1*2] = 1.0;
      t[1+1*2] = 0.0;
      t[0+2*2] = 0.5;
      t[1+2*2] = 10.0;
    }
    else if ( i == 4 )
    {
      t[0+0*2] = 0.0;
      t[1+0*2] = 0.0;
      t[0+1*2] = 1.0;
      t[1+1*2] = 0.0;
      t[0+2*2] = 10.0;
      t[1+2*2] = 2.0;
    }

    r8mat_transpose_print ( DIM_NUM, 3, t, "  The triangle" );

    center = triangle_circumcenter_2d ( t );

    r8vec_print ( DIM_NUM, center, "  Circumcenter" );

    delete [] center;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests TRIANGULATION_PLOT_EPS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define NODE_NUM 9
# define DIM_NUM 2
# define TRI_NUM 12

  char file_name[80] = "triangulation_plot.eps";
  double g_xy[DIM_NUM*NODE_NUM] = {
       0.0, 0.0,
       0.0, 1.0,
       0.2, 0.5,
       0.3, 0.6,
       0.4, 0.5,
       0.6, 0.4,
       0.6, 0.5,
       1.0, 0.0,
       1.0, 1.0 };
  int nod_tri[TRI_NUM*3] = {
       2, 1, 3,
       3, 1, 6,
       2, 3, 4,
       4, 3, 5,
       7, 4, 5,
       5, 3, 6,
       7, 5, 6,
       9, 4, 7,
       6, 1, 8,
       7, 6, 8,
       7, 8, 9,
       2, 4, 9 };

  cout << "\n";
  cout << "TEST07\n";
  cout << "  TRIANGULATION_PLOT_EPS can plot a triangulation.\n";

  triangulation_plot_eps ( file_name, NODE_NUM, g_xy, TRI_NUM, nod_tri );

  cout << "\n";
  cout << "  TRIANGULATION_PLOT_EPS has created an\n";
  cout << "  Encapsulated PostScript file (EPS) containing\n";
  cout << "  an image of the triangulation.\n";
  cout << "\n";
  cout << "  This file is called " << file_name << "\n";

  return;
# undef NODE_NUM
# undef DIM_NUM
# undef TRI_NUM
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests TRIANGULATION_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define NODE_NUM 9
# define DIM_NUM 2
# define TRI_NUM 12

  double g_xy[DIM_NUM*NODE_NUM] = {
       0.0, 0.0,
       0.0, 1.0,
       0.2, 0.5,
       0.3, 0.6,
       0.4, 0.5,
       0.6, 0.4,
       0.6, 0.5,
       1.0, 0.0,
       1.0, 1.0 };
  int nod_tri[TRI_NUM*3] = {
       2, 1, 3,
       3, 1, 6,
       2, 3, 4,
       4, 3, 5,
       7, 4, 5,
       5, 3, 6,
       7, 5, 6,
       9, 4, 7,
       6, 1, 8,
       7, 6, 8,
       7, 8, 9,
       2, 4, 9 };
  int triangle_neighbor[TRI_NUM*3] = {
       -28,   2,  3,
         1,   9,  6,
         1,   4, 12,
         3,   6,  5,
         8,   4,  7,
         4,   2,  7,
         5,   6, 10,
        12,   5, 11,
         2, -34, 10,
         7,   9, 11,
        10, -38,  8,
         3,   8, -3 };

  cout << "\n";
  cout << "TEST08\n";
  cout << "  TRIANGULATION_PRINT prints out a triangulation.\n";

  triangulation_print ( NODE_NUM, g_xy, TRI_NUM, nod_tri, triangle_neighbor );

  return;
# undef NODE_NUM
# undef DIM_NUM
# undef TRI_NUM
}
