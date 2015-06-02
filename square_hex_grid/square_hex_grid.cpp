# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <fstream>
# include <cstring>

using namespace std;

# include "square_hex_grid.hpp"

//****************************************************************************80

void box_print_2d ( double box[] )

//****************************************************************************80
//
//  Purpose:
//
//    BOX_PRINT_2D prints information about a coordinate box in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double BOX[2*2], the coordinates of the lower left
//    and upper right corners of the box.
//
{
  cout << "\n";
  cout << "  Coordinate box:\n";
  cout << "  "
       << setw(12) << box[0+0*2] << " <= X <= "
       << setw(12) << box[0+1*2] << "\n";
  cout << "  "
       << setw(12) << box[1+0*2] << " <= Y <= "
       << setw(12) << box[1+1*2] << "\n";

  return;
}
//****************************************************************************80

char digit_to_ch ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_TO_CH returns the base 10 digit character corresponding to a digit.
//
//  Example:
//
//     I     C
//   -----  ---
//     0    '0'
//     1    '1'
//   ...    ...
//     9    '9'
//    10    '*'
//   -83    '*'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the digit, which should be between 0 and 9.
//
//    Output, char DIGIT_TO_CH, the appropriate character '0' through '9' or '*'.
//
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else
  {
    c = '*';
  }

  return c;
}
//****************************************************************************80

void hex_grid_01_approximate_h ( double h_goal, int *nodes_per_layer,
  double *h )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_01_APPROXIMATE_H seeks a unit square hex grid with spacing H.
//
//  Discussion:
//
//    The parameter NODES_PER_LAYER controls the number of nodes and the
//    grid spacing, but in a somewhat obscure way.  This routine experiments
//    with various values until it is convinced it has the value
//    of NODES_PER_LAYER that produces a grid spacing that is no
//    no greater than H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double H_GOAL, the desired grid spacing.
//
//    Output, int *NODES_PER_LAYER, the number of nodes per layer
//    which produces a mesh with grid spacing H_GOAL or less.
//
//    Output, double *H, the actual grid spacing.
//
{
  int nodes_per_layer2;

  if ( h_goal <= 0.0 )
  {
    cout << "\n";
    cout << "HEX_GRID_01_APPROXIMATE_H - Fatal error!\n";
    cout << "  Illegal input value of H_GOAL = " << h_goal << "\n";
    exit ( 1 );
  }

  *nodes_per_layer = 1 + ( int ) ( 1.0 / h_goal );
//
//  Check whether roundoff means we could use one less node per layer.
//
  if ( 2 < *nodes_per_layer )
  {
    nodes_per_layer2 = *nodes_per_layer - 1;
    *h = 1.0 / ( double ) ( nodes_per_layer2 - 1 );

    if ( *h <= h_goal )
    {
      *nodes_per_layer = nodes_per_layer2;
      return;
    }

  }

  *h = 1.0 / ( double ) ( *nodes_per_layer - 1 );

  return;
}
//****************************************************************************80

void hex_grid_01_approximate_n ( int n_goal, int *nodes_per_layer, int *n )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_01_APPROXIMATE_N seeks a unit square hex grid of about N nodes.
//
//  Discussion:
//
//    The parameter NODES_PER_LAYER controls the number of nodes, but
//    in a somewhat obscure way.  This routine experiments with various
//    values until it is convinced it has the value of NODES_PER_LAYER
//    that comes as close as possible to producing N nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N_GOAL, the desired number of nodes.
//
//    Output, int *NODES_PER_LAYER, the number of nodes per layer
//    which produces a mesh with about N_GOAL nodes.
//
//    Output, int *N, the number of nodes in the mesh.
//
{
  int n_high;
  int n_low;
  int nodes_per_layer_high;
  int nodes_per_layer_low;

  if ( n_goal <= 1 )
  {
    cout << "\n";
    cout << "HEX_GRID_01_APPROXIMATE_N - Fatal error!\n";
    cout << "  Illegal input value of N_GOAL = " << n_goal << "\n";
    exit ( 1 );
  }

  nodes_per_layer_low = 0;
  n_low = 0;

  *nodes_per_layer = ( int ) ( 0.5 + sqrt ( ( double ) ( n_goal ) ) );

  nodes_per_layer_high = n_goal;
  n_high = n_goal * n_goal;

  for ( ; ; )
  {
    *n = hex_grid_01_n ( *nodes_per_layer );

    if ( *n == n_goal )
    {
      break;
    }

    if ( *n < n_goal )
    {
      nodes_per_layer_low = *nodes_per_layer;
      n_low = *n;
    }
    else
    {
      nodes_per_layer_high = *nodes_per_layer;
      n_high = *n;
    }

    if ( nodes_per_layer_low + 1 == nodes_per_layer_high )
    {
      if ( *n - n_low <= n_high - *n )
      {
        *nodes_per_layer = nodes_per_layer_high;
        *n = n_high;
      }
      else
      {
        *nodes_per_layer = nodes_per_layer_low;
        *n = n_low;
      }
      break;
    }

    *nodes_per_layer = ( nodes_per_layer_low + nodes_per_layer_high ) / 2;

  }

  return;
}
//****************************************************************************80

void hex_grid_01_h ( int nodes_per_layer, double *hx, double *hy )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_01_H computes the unit square hex grid spacings.
//
//  Discussion:
//
//    This routine determines the values of HX and HY from
//    the fundamental hexagonal grid parameter NODES_PER_LAYER.
//
//    A hexagonal grid is defined in the unit square [0,1] x [0,1].
//
//    All nodes of the grid lie on one of LAYERS horizontal lines.
//    The first of these lines is the X axis, and each successive
//    line is HY units higher.
//
//    On all the odd numbered lines, there are NODES_PER_LAYER points,
//    equally spaced from 0 to 1, with a spacing of HX.
//
//    On the even numbered lines, there are NODES_PER_LAYER-1 points,
//    whose values are the midpoints of successive intervals on
//    an odd numbered line.  (The grid is staggered).
//
//    HY = HX * sqrt ( 3 ) / 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODES_PER_LAYER, the number of grid points on the first
//    horizontal layer of points.
//
//    Output, double *HX, the spacing between grid points
//    on a horizontal line.
//
//    Output, double *HY, the spacing between horizontal lines.
//
{
  if ( nodes_per_layer < 1 )
  {
    *hx = 0.0;
    *hy = 0.0;
  }
  else if ( nodes_per_layer == 1 )
  {
    *hx = 1.0;
    *hy = 1.0;
  }
  else
  {
    *hx = 1.0 / ( double ) ( nodes_per_layer - 1 );
    *hy = ( *hx ) * sqrt ( 3.0 ) / 2.0;
  }

  return;
}
//****************************************************************************80

int hex_grid_01_layers ( int nodes_per_layer )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_01_LAYERS computes the unit square hex grid column width.
//
//  Discussion:
//
//    This routine determines the value of LAYERS, the number of
//    layers, from the fundamental hexagonal grid parameter NODES_PER_LAYER.
//
//    A hexagonal grid is defined in the unit square [0,1] x [0,1].
//
//    All nodes of the grid lie on one of LAYERS horizontal lines.
//    The first of these lines is the X axis, and each successive
//    line is HY units higher.
//
//    On all the odd numbered lines, there are NODES_PER_LAYER points,
//    equally spaced from 0 to 1, with a spacing of HX.
//
//    On the even numbered lines, there are NODES_PER_LAYER-1 points,
//    whose values are the midpoints of successive intervals on
//    an odd numbered line.  (The grid is staggered).
//
//    HY = HX * sqrt ( 3 ) / 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODES_PER_LAYER, the number of grid points on the first
//    horizontal layer of points.
//
//    Output, int HEX_GRID_01_LAYERS, the number of horizontal layers.
//
{
  double hx;
  double hy;
  int layers;

  if ( nodes_per_layer < 1 )
  {
    layers = 0;
  }
  else if ( nodes_per_layer == 1 )
  {
    layers = 1;
  }
  else
  {
    hx = 1.0 / ( double ) ( nodes_per_layer - 1 );
    hy = sqrt ( 3.0 ) * hx / 2.0;
    layers = 1 + ( int ) ( 1.0 / hy );
  }

  return layers;
}
//****************************************************************************80

int hex_grid_01_n ( int nodes_per_layer )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_01_N computes the number of unit square hex grid points.
//
//  Discussion:
//
//    This routine determines the value of N, the number of
//    hex grid points, from the fundamental hexagonal grid
//    parameter NODES_PER_LAYER.
//
//    A hexagonal grid is defined in the unit square [0,1] x [0,1].
//
//    All nodes of the grid lie on one of LAYERS horizontal lines.
//    The first of these lines is the X axis, and each successive
//    line is HY units higher.
//
//    On all the odd numbered lines, there are NODES_PER_LAYER points,
//    equally spaced from 0 to 1, with a spacing of HX.
//
//    On the even numbered lines, there are NODES_PER_LAYER-1 points,
//    whose values are the midpoints of successive intervals on
//    an odd numbered line.  (The grid is staggered).
//
//    HY = HX * sqrt ( 3 ) / 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODES_PER_LAYER, the number of grid points on the first
//    horizontal layer of points.
//
//    Output, int HEX_GRID_01_N, the number of hex grid points.
//
{
  int n;
  int layers;

  if ( nodes_per_layer < 1 )
  {
    n = 0;
  }
  else if ( nodes_per_layer == 1 )
  {
    n = 1;
  }
  else
  {
    layers = hex_grid_01_layers ( nodes_per_layer );

    n = nodes_per_layer       * ( ( layers + 1 ) / 2 ) +
      ( nodes_per_layer - 1 ) * ( ( layers     ) / 2 );
  }

  return n;
}
//****************************************************************************80

double *hex_grid_01_points ( int nodes_per_layer, int layers, int n )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_01_POINTS returns unit square hex grid points.
//
//  Discussion:
//
//    This routine determines the coordinates of the elements of
//    a hexagonal grid in the unit square.
//
//    A hexagonal grid is defined in the unit square [0,1] x [0,1].
//
//    All nodes of the grid lie on one of LAYERS horizontal lines.
//    The first of these lines is the X axis, and each successive
//    line is HY units higher.
//
//    On all the odd numbered lines, there are NODES_PER_LAYER points,
//    equally spaced from 0 to 1, with a spacing of HX.
//
//    On the even numbered lines, there are NODES_PER_LAYER-1 points,
//    whose values are the midpoints of successive intervals on
//    an odd numbered line.  (The grid is staggered).
//
//    HY = HX * sqrt ( 3 ) / 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODES_PER_LAYER, the number of grid points on the first
//    horizontal layer of points.
//
//    Input, int LAYERS, the number of horizontal layers.
//
//    Input, int N, the total number of hex grid points.
//
//    Output, float P[2*N], the coordinates of the
//    mesh points, listed one horizontal layer at a time.
//
{
  double hx;
  double hy;
  int i;
  int j;
  int jmod;
  int k;
  int ndim = 2;
  double *p;
  double x;
  double y;

  if ( nodes_per_layer < 1 )
  {
    return NULL;
  }

  p = new double[ndim*n];

  if ( nodes_per_layer == 1 )
  {
    for ( i = 0; i < ndim; i++ )
    {
      p[i+0*2] = 0.5;
    }
    return p;
  }

  hex_grid_01_h ( nodes_per_layer, &hx, &hy );

  k = 0;

  for  ( j = 1; j <= layers; j++ )
  {
    y = hy * ( double ) ( j - 1 );

    jmod = j % 2;

    if ( jmod == 1 )
    {
      for ( i = 1; i <= nodes_per_layer; i++ )
      {
        x = ( double ) ( i - 1 ) / ( double ) ( nodes_per_layer - 1 );
        k = k + 1;
        if ( k <= n )
        {
          p[0+(k-1)*2] = x;
          p[1+(k-1)*2] = y;
        }
      }
    }
    else
    {
      for ( i = 1; i <= nodes_per_layer-1; i++ )
      {
        x = ( double ) ( 2 * i - 1 ) / ( double ) ( 2 * nodes_per_layer - 2 );
        k = k + 1;
        if ( k <= n )
        {
          p[0+(k-1)*2] = x;
          p[1+(k-1)*2] = y;
        }
      }

    }

  }

  return p;
}
//****************************************************************************80

void hex_grid_approximate_h ( double box[], double h_goal, int *nodes_per_layer,
  double *h )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_APPROXIMATE_H seeks a hex grid with spacing H.
//
//  Discussion:
//
//    The parameter NODES_PER_LAYER controls the number of nodes and the
//    grid spacing, but in a somewhat obscure way.  This routine experiments
//    with various values until it is convinced it has the value
//    of NODES_PER_LAYER that produces a grid spacing that is no
//    no greater than H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double BOX[2*2], the lower and upper corners
//    of the rectangular region.
//
//    Input, double H_GOAL, the desired grid spacing.
//
//    Output, int *NODES_PER_LAYER, the number of nodes per layer
//    which produces a mesh with grid spacing H_GOAL or less.
//
//    Output, double *H, the actual grid spacing.
//
{
  int nodes_per_layer2;

  if ( h_goal <= 0.0 )
  {
    cout << "\n";
    cout << "HEX_GRID_APPROXIMATE_H - Fatal error!\n";
    cout << "  Illegal input value of H_GOAL = " << h_goal <<"\n";
    exit ( 1 );
  }

  *nodes_per_layer = 1 + ( int ) ( ( box[0+1*2] - box[0+0*2] ) / h_goal );
//
//  Check whether roundoff means we could use one less node per layer.
//
  if ( 2 < *nodes_per_layer )
  {
    nodes_per_layer2 = *nodes_per_layer - 1;
    *h = ( box[0+1*2] - box[0+0*2] ) / ( double ) ( nodes_per_layer2 - 1 );

    if ( *h <= h_goal )
    {
      *nodes_per_layer = nodes_per_layer2;
      return;
    }

  }

  *h = ( box[0+1*2] - box[0+0*2] ) / ( double ) ( *nodes_per_layer - 1 );

  return;
}
//****************************************************************************80

void hex_grid_approximate_n ( double box[], int n_goal, int *nodes_per_layer,
  int *n )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_APPROXIMATE_N seeks a hex grid of about N nodes.
//
//  Discussion:
//
//    The parameter NODES_PER_LAYER controls the number of nodes, but
//    in a somewhat obscure way.  This routine experiments with various
//    values until it is convinced it has the value of NODES_PER_LAYER
//    that comes as close as possible to producing N nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double BOX[2*2], the lower and upper corners
//    of the rectangular region.
//
//    Input, int N_GOAL, the desired number of nodes.
//
//    Output, int *NODES_PER_LAYER, the number of nodes per layer
//    which produces a mesh with about N_GOAL nodes.
//
//    Output, int *N, the number of nodes in the mesh.
//
{
  int n_high;
  int n_low;
  int nodes_per_layer_high;
  int nodes_per_layer_low;

  if ( n_goal <= 1 )
  {
    cout << "\n";
    cout << "HEX_GRID_APPROXIMATE_N - Fatal error!\n";
    cout << "  Illegal input value of N_GOAL = " << n_goal << "\n";
    exit ( 1 );
  }

  nodes_per_layer_low = 0;
  n_low = 0;

  *nodes_per_layer = ( int ) ( 0.5 + sqrt ( ( double ) ( n_goal ) ) );

  nodes_per_layer_high = n_goal;
  n_high = n_goal * n_goal;

  for ( ; ; )
  {
    *n = hex_grid_n ( *nodes_per_layer, box );

    if ( *n == n_goal )
    {
      break;
    }

    if ( *n < n_goal )
    {
      nodes_per_layer_low = *nodes_per_layer;
      n_low = *n;
    }
    else
    {
      nodes_per_layer_high = *nodes_per_layer;
      n_high = *n;
    }

    if ( nodes_per_layer_low + 1 == nodes_per_layer_high )
    {
      if ( *n - n_low <= n_high - *n )
      {
        *nodes_per_layer = nodes_per_layer_high;
        *n = n_high;
      }
      else
      {
        *nodes_per_layer = nodes_per_layer_low;
        *n = n_low;
      }
      break;
    }

    *nodes_per_layer = ( nodes_per_layer_low + nodes_per_layer_high ) / 2;

  }

  return;
}
//****************************************************************************80

void hex_grid_h ( int nodes_per_layer, double box[], double *hx, double *hy )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_H computes the coordinate box hex grid spacings.
//
//  Discussion:
//
//    This routine determines the values of HX and HY from
//    the fundamental hexagonal grid parameter NODES_PER_LAYER.
//
//    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
//
//    All nodes of the grid lie on one of LAYERS horizontal lines.
//    The first of these lines is from (A,C) to (B,C), and each
//    successive line is HY units higher.
//
//    On all the odd numbered lines, there are NODES_PER_LAYER points,
//    equally spaced from A to B, with a spacing of HX.
//
//    On the even numbered lines, there are NODES_PER_LAYER-1 points,
//    whose values are the midpoints of successive intervals on
//    an odd numbered line.  (The grid is staggered).
//
//    HY = HX * sqrt ( 3 ) / 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODES_PER_LAYER, the number of grid points on the first
//    horizontal layer of points.
//
//    Input, double BOX[2*2], the values of A, B, C and D
//    that define the coordinate box.
//
//    Output, double *HX, the spacing between grid points
//    on a horizontal line.
//
//    Output, double *HY, the spacing between horizontal lines.
//
{
  if ( nodes_per_layer < 1 )
  {
    *hx = 0.0;
    *hy = 0.0;
  }
  else if ( nodes_per_layer == 1 )
  {
    *hx = box[0+1*2] - box[0+0*2];
    *hy = box[1+1*2] - box[1+0*2];
  }
  else
  {
    *hx = ( box[0+1*2] - box[0+0*2] ) / ( double ) ( nodes_per_layer - 1 );
    *hy = ( *hx ) * sqrt ( 3.0 ) / 2.0;
  }

  return;
}
//****************************************************************************80

int hex_grid_layers ( int nodes_per_layer, double box[] )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_LAYERS computes the coordinate box hex grid column width.
//
//  Discussion:
//
//    This routine determines the value of LAYERS, the number of
//    layers, from the fundamental hexagonal grid parameter NODES_PER_LAYER.
//
//    A hexagonal grid is defined in a coordinate box [A,B] x [C,D].
//
//    All nodes of the grid lie on one of LAYERS horizontal lines.
//    The first of these lines is from (A,C) to (B,C), and each
//    successive line is HY units higher.
//
//    On all the odd numbered lines, there are NODES_PER_LAYER points,
//    equally spaced from A to B, with a spacing of HX.
//
//    On the even numbered lines, there are NODES_PER_LAYER-1 points,
//    whose values are the midpoints of successive intervals on
//    an odd numbered line.  (The grid is staggered).
//
//    HY = HX * sqrt ( 3 ) / 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODES_PER_LAYER, the number of grid points on the first
//    horizontal layer of points.
//
//    Input, double BOX[2*2], the values of A, B, C and D
//    that define the coordinate box.
//
//    Output, int HEX_GRID_LAYERS, the number of horizontal layers.
//
{
  double hx;
  double hy;
  int layers;

  if ( nodes_per_layer < 1 )
  {
    layers = 0;
  }
  else if ( nodes_per_layer == 1 )
  {
    layers = 1;
  }
  else
  {
    hx = ( box[0+1*2] - box[0+0*2] ) / ( double ) ( nodes_per_layer - 1 );
    hy = sqrt ( 3.0 ) * hx / 2.0;
    layers = 1 + ( int ) ( ( box[1+1*2] - box[1+0*2] ) / hy );
  }

  return layers;
}
//****************************************************************************80

int hex_grid_n ( int nodes_per_layer, double box[] )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_N computes the number of coordinate box hex grid points.
//
//  Discussion:
//
//    This routine determines the value of N, the number of
//    hex grid points, from the fundamental hexagonal grid
//    parameter NODES_PER_LAYER.
//
//    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
//
//    All nodes of the grid lie on one of LAYERS horizontal lines.
//    The first of these lines is from (A,C) to (B,C), and each
//    successive line is HY units higher.
//
//    On all the odd numbered lines, there are NODES_PER_LAYER points,
//    equally spaced from A to B, with a spacing of HX.
//
//    On the even numbered lines, there are NODES_PER_LAYER-1 points,
//    whose values are the midpoints of successive intervals on
//    an odd numbered line.  (The grid is staggered).
//
//    HY = HX * sqrt ( 3 ) / 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODES_PER_LAYER, the number of grid points on the first
//    horizontal layer of points.
//
//    Input, double BOX[2*2], the values of A, B, C and D
//    that define the coordinate box.
//
//    Output, int HEX_GRID_N, the number of hex grid points.
//
{
  int n;
  int layers;

  if ( nodes_per_layer < 1 )
  {
    n = 0;
  }
  else if ( nodes_per_layer == 1 )
  {
    n = 1;
  }
  else
  {
    layers = hex_grid_layers ( nodes_per_layer, box );

    n = nodes_per_layer       * ( ( layers + 1 ) / 2 ) +
      ( nodes_per_layer - 1 ) * ( ( layers     ) / 2 );
  }

  return n;
}
//****************************************************************************80

double *hex_grid_points ( int nodes_per_layer, int layers, int n, double box[] )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_POINTS returns coordinate box hex grid points.
//
//  Discussion:
//
//    This routine determines the coordinates of the elements of
//    a hexagonal grid in the unit square.
//
//    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
//
//    All nodes of the grid lie on one of LAYERS horizontal lines.
//    The first of these lines is from (A,C) to (B,C), and each
//    successive line is HY units higher.
//
//    On all the odd numbered lines, there are NODES_PER_LAYER points,
//    equally spaced from A to B, with a spacing of HX.
//
//    On the even numbered lines, there are NODES_PER_LAYER-1 points,
//    whose values are the midpoints of successive intervals on
//    an odd numbered line.  (The grid is staggered).
//
//    HY = HX * sqrt ( 3 ) / 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODES_PER_LAYER, the number of grid points on the first
//    horizontal layer of points.
//
//    Input, int LAYERS, the number of horizontal layers.
//
//    Input, int N, the total number of hex grid points.
//
//    Input, double BOX[2*2], the values of A, B, C and D
//    that define the coordinate box.
//
//    Output, double HEX_GRID_POINTS[2*N], the coordinates of the
//    mesh points, listed one horizontal layer at a time.
//
{
  double hx;
  double hy;
  int i;
  int j;
  int jmod;
  int k;
  int ndim = 2;
  double *p;
  double x;
  double y;

  if ( nodes_per_layer < 1 )
  {
    return NULL;
  }

  p = new double[n*ndim];

  if ( nodes_per_layer == 1 )
  {
    for ( i = 0; i < 2; i++ )
    {
      p[i+0*2] = ( box[i+0*2] + box[i+1*2] ) / 2.0;
    }
    return p;
  }

  hex_grid_h ( nodes_per_layer, box, &hx, &hy );

  k = 0;

  for ( j = 1; j <= layers; j++ )
  {
    y = box[1+0*2] + hy * ( double ) ( j - 1 );

    jmod = j % 2;

    if ( jmod == 1 )
    {
      for ( i = 1; i <= nodes_per_layer; i++ )
      {
        x = box[0+0*2] + ( box[0+1*2] - box[0+0*2] ) * ( double ) ( i - 1 )
          / ( double ) ( nodes_per_layer - 1 );
        k = k + 1;
        if ( k <= n )
        {
          p[0+(k-1)*2] = x;
          p[1+(k-1)*2] = y;
        }
      }
    }
    else
    {
      for ( i = 1; i <= nodes_per_layer-1; i++ )
      {
        x = box[0+0*2] + ( box[0+1*2] - box[0+0*2] ) * ( double ) ( 2 * i - 1 )
          / ( double ) ( 2 * nodes_per_layer - 2 );
        k = k + 1;
        if ( k <= n )
        {
          p[0+(k-1)*2] = x;
          p[1+(k-1)*2] = y;
        }
      }

    }

  }

  return p;
}
//****************************************************************************80

void hex_grid_write ( int n, int nodes_per_layer, int layers, double hx,
  double hy, double box[], double r[], char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_GRID_WRITE writes a hex grid dataset to a file.
//
//  Discussion:
//
//    The initial lines of the file are comments, which begin with a
//    '#' character.
//
//    Thereafter, each line of the file contains the M-dimensional
//    components of the next entry of the sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, int NODES_PER_LAYER, the number of nodes in the first layer.
//
//    Input, int LAYERS, the number of layers.
//
//    Input, double HX, HY, the row and column spacings.
//
//    Input, double BOX[2*2], the values of A, B, C and D
//    that define the coordinate box.
//
//    Input, double R[2*N], the points.
//
//    Input, char *FILE_OUT_NAME, the name of
//    the output file.
//
{
  ofstream file_out;
  int i;
  int j;
  char *s;

  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "HEX_GRID_WRITE - Fatal error!\n";
    cout << "  Could not open the output file.\n";
    exit ( 1 );
  }

  s = timestring ( );

  file_out << "#  " << file_out_name << "\n";
  file_out << "#  created by routine HEX_GRID_WRITE.C" << "\n";
  file_out << "#  at " << s << "\n";
  file_out << "#\n";
  file_out << "#  Spatial dimension M =     "
           << setw(12) <<  2 << "\n";
  file_out << "#  Number of points N =      "
           << setw(12) <<  n << "\n";
  file_out << "#  NODES_PER_LAYER =         "
           << setw(12) <<  nodes_per_layer << "\n";
  file_out << "#  Number of rows LAYERS =   "
           << setw(12) <<  layers << "\n";
  file_out << "#  Coordinate box X(1) =     "
           << setw(12) <<  box[0+0*2] << "\n";
  file_out << "#  Coordinate box X(2) =     "
           << setw(12) <<  box[0+1*2] << "\n";
  file_out << "#  Coordinate box Y(1) =     "
           << setw(12) <<  box[1+0*2] << "\n";
  file_out << "#  Coordinate box Y(2) =     "
           << setw(12) <<  box[1+1*2] << "\n";
  file_out << "#  Node spacing HX =         "
           << setw(12) <<  hx << "\n";
  file_out << "#  Layer spacing HY =        "
           << setw(12) <<  hy << "\n";
  file_out << "#  EPSILON (unit roundoff) = "
           << setw(12) <<  r8_epsilon ( ) << "\n";

  file_out << "#\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      file_out << setw(10) << r[i+j*2] << "  ";
    }
    file_out << "\n";
  }

  file_out.close ( );

  delete [] s;

  return;
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
//    It should be the case that 10^I4_LOG_10(I) <= |I| < 10^(I4_LOG_10(I)+1).
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

char *i4_to_s ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_S converts an I4 to a string.
//
//  Example:
//
//    INTVAL  S
//
//         1  1
//        -1  -1
//         0  0
//      1952  1952
//    123456  123456
//   1234567  1234567
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, an integer to be converted.
//
//    Output, char *I4_TO_S, the representation of the integer.
//
{
  int digit;
  int j;
  int length;
  int ten_power;
  char *s;

  length = i4_log_10 ( i );

  ten_power = ( int ) pow ( ( double ) 10, ( double ) length );

  if ( i < 0 )
  {
    length = length + 1;
  }
//
//  Add one position for the trailing null.
//
  length = length + 1;

  s = new char[length];

  if ( i == 0 )
  {
    s[0] = '0';
    s[1] = '\0';
    return s;
  }
//
//  Now take care of the sign.
//
  j = 0;
  if ( i < 0 )
  {
    s[j] = '-';
    j = j + 1;
    i = abs ( i );
  }
//
//  Find the leading digit of I, strip it off, and stick it into the string.
//
  while ( 0 < ten_power )
  {
    digit = i / ten_power;
    s[j] = digit_to_ch ( digit );
    j = j + 1;
    i = i - digit * ten_power;
    ten_power = ten_power / 10;
  }
//
//  Tack on the trailing NULL.
//
  s[j] = '\0';
  j = j + 1;

  return s;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
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
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

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
//    May 31 2001 09:45:54 AM
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

char *timestring ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTRING returns the current YMDHMS date as a string.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
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
