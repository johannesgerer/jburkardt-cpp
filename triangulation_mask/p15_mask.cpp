//****************************************************************************80

bool triangle_mask ( int dim_num, int triangle_order, int nodes[], 
  double coord[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_MASK is a user routine which masks triangles.
//
//  Discussion:
//
//    The region to be considered is the union of two rectangles.
//    The first is  -8 <= X <= 2, -1 <= Y <= 0,
//    the second is -2 <= X <= 8,  0 <= Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int TRIANGLE_ORDER, the number of nodes in the triangle.
//
//    Input, int NODES[TRIANGLE_ORDER], the indices of the nodes.
//
//    Input, double COORD[DIM_NUM*TRIANGLE_ORDER], the coordinates
//    of the nodes.
//
//    Output, bool TRIANGLE_MASK, is TRUE if the triangle should be discarded,
//    and FALSE if the triangle should be retained.
//
{
  double centroid[2];
  int dim;
  bool mask;
  int order;
//
//  Compute the centroid.
//
  for ( dim = 0; dim < 2; dim++ )
  {
    centroid[dim] = 0.0;
    for ( order = 0; order < triangle_order; order++ )
    {
      centroid[dim] = centroid[dim] + coord[dim+order*dim_num];
    }
    centroid[dim] = centroid[dim] / ( double ) ( triangle_order );
  }
//
//  MASK = The centroid is outside the region.
//
  if (      -8.0 <= centroid[0]        &&
                    centroid[0] <= 2.0 &&
            -1.0 <= centroid[1]        &&
                    centroid[1] <= 0.0 )
  {
    mask = false;
  }
  else if ( -2.0 <= centroid[0]        &&
                    centroid[0] <= 8.0 &&
             0.0 <= centroid[1]        &&
                    centroid[1] <= 1.0 )
  {
    mask = false;
  }
  else
  {
    mask = true;
  }

  return mask;
}
