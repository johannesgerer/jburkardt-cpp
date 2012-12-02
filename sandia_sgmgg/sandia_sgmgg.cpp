# include "sandia_rules.hpp"
# include "sandia_sgmgg.hpp"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <ctime>

namespace webbur
{
//****************************************************************************80

void sandia_sgmgg_coef_inc2 ( int m, int n1, int s1[], int c1[],
  int s2[], int c3[] )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGG_COEF_INC2 computes tentative coefficient changes.
//
//  Discussion:
//
//    An active set S1 of N1 sparse grid indices is given, each of
//    size M.
//
//    The coefficient C1 of each sparse grid index is also given.
//
//    A candidate sparse grid index S2 is provided.
//
//    This function determines the N+1 coefficients that would be
//    appropriate if the candidate S2 was added to the active set
//    as the (N+1)-st item.
//
//    During the calculation, we may try to update coefficients of inactive 
//    index sets.  By the end of the calculation, all these inactive index
//    sets should have accumulated total coefficients of 0 again.  As a check,
//    we temporarily set aside space for these objects, and check, at the end,
//    that the coefficients are zero.
//
//  Example:
//
//    Input:
//
//      +1 * {0,2}
//      -1 * {0,1}  +1 * {1,1}
//                  -1 * {1,0}  +1 * {2,0}
//
//    Add {3,0}
//
//    Output:
//
//      +1 * {0,2}
//      -1 * {0,1}  +1 * {1,1}
//                  -1 * {1,0}   0 * {2,0}  +1 * {3,0}
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the vector.
//
//    Input, int N1, the number of points in the active set.
//
//    Input, int S1[M*N1], the indices for the active set.
//
//    Input, int C1[N1], the coefficients for the active set.
//
//    Input, int S2[M], the indices for the candidate.
//
//    Output, int C3[N1+1], the coefficients for the active set
//    plus the candidate.
//
{
  int *c4;
  int i;
  int i2;
  int j;
  int j2;
  int k;
  int n4;
  int *s;
  int *s4;
//
//  Initialize the inactive data.
//
  n4 = 0;

  c4 = new int[n1];
  for ( i = 0; i < n1; i++ )
  {
    c4[i] = 0;
  }

  s4 = new int[m*n1];
  for ( j = 0; j < n1; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      s4[i+j*m] = 0;
    }
  }
//
//  Copy C1.
//
  for ( j = 0; j < n1; j++ )
  {
    c3[j] = c1[j];
  }
  c3[n1] = 1;
//
//  Consider the effect of the new item S2 on each of the current
//  items in the active set S1.
//
  s = new int[m];

  for ( j = 0; j < n1; j++ )
  {
//
//  Determine S, the element-wise minimum of the J-th item in S1 versus S2.
//
    k = j;

    for ( i = 0; i < m; i++ )
    {
      if ( s2[i] < s1[i+j*m] )
      {
        s[i] = s2[i];
        k = -1;
      }
      else
      {
        s[i] = s1[i+j*m];
      }
    }
//
//  If S = S1[*,J], then K is J.
//
    if ( k != -1 )
    {
      c3[k] = c3[k] - c1[j];
    }
//
//  If S is equal to an element of the active set, we set K to that index.
//
    else
    {
      for ( j2 = 0; j2 < n1; j2++ )
      {
        k = j2;
        for ( i2 = 0; i2 < m; i2++ )
        {
          if ( s1[i2+j2*m] != s[i2] )
          {
             k = -1;
             break;
          }
        }
        if ( k != -1 )
        {
           c3[k] = c3[k] - c1[j];
           break;
        }
      }
    }
//
//  If S is equal to an element of the inactive set, set K to that index.
//
    if ( k == -1 )
    {
      for ( j2 = 0; j2 < n4; j2++ )
      {
        k = j2;

        for ( i2 = 0; i2 < m; i2++ )
        {
          if ( s4[i2+j2*m] != s[i2] )
          {
            k = - 1;
            break;
          }
        }

        if ( k != - 1 )
        {
          c4[k] = c4[k] - c1[j];
          break;
        }
      }
    }
//
//  S is not equal to S1(*,J), or any element of S1, or any element of S4.
//  Add S to the set of elements S4.
//
    if ( k == -1 )
    {
      k = n4;
      c4[k] = 0;
      for ( i = 0; i < m; i++ )
      {
        s4[i+k*m] = s[i];
      }
      c4[k] = c4[k] - c1[j];
      n4 = n4 + 1;
    }
  }
//
//  At the end, the C4(1:N4) should all be zero.
//
  for ( i = 0; i < n4; i++ )
  {
    if ( c4[k] != 0 )
    {
      std::cerr << "\n";
      std::cerr << "SANDIA_SGMGG_COEF_INC2 - Fatal error!\n";
      std::cerr << "  Some inactive indices were assigned a nonzero coefficient.\n";
      webbur::i4mat_transpose_print ( m, n4, s4, "  S4:" );
      webbur::i4vec_print ( n4, c4, "  C4:" );
      std::exit ( 1 );
    }
  }

  delete [] c4;
  delete [] s;
  delete [] s4;

  return;
}
//****************************************************************************80

void sandia_sgmgg_coef_naive ( int dim_num, int point_num, int sparse_index[],
  int coef[] )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGG_COEF_NAIVE returns the combinatorial coefficients.
//
//  Discussion:
//
//    The coefficient of point I is calculated as follows:
//
//    *) point J is a "neighbor" of point I if every entry of the sparse
//       index for point J is either equal to, or 1 greater than, the
//       corresponding entry of the sparse index of point I.
//
//    *) If point J is a neighbor of point I, then it contributes
//       (-1)^D to the coefficient, where D is the sum of the differences
//       between the sparse indices of point I and point J.
//
//    This is a completely naive implementation of the calculation,
//    intended simply as a demonstration for small examples.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the vector.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int SPARSE_INDEX[DIM_NUM*POINT_NUM],
//    the indices that define the points.
//
//    Output, int COEF[POINT_NUM], the coefficients.
//
{
  int dif;
  int i;
  int j;
  int j1;
  int j2;
  bool neighbor;
  int term;

  for ( j = 0; j < point_num; j++ )
  {
    coef[j] = 0;
  }
  for ( j1 = 0; j1 < point_num; j1++ )
  {
    for ( j2 = 0; j2 < point_num; j2++ )
    {
      neighbor = true;
      term = + 1;
      for ( i = 0; i < dim_num; i++ )
      {
        dif = sparse_index[i+j2*dim_num] - sparse_index[i+j1*dim_num];

        if ( dif == 0 )
        {
        }
        else if ( dif == 1 )
        {
          term = - term;
        }
        else
        {
          neighbor = false;
          break;
        }
      }
      if ( neighbor )
      {
        coef[j1] = coef[j1] + term;
      }
    }
  }
  return;
}
//****************************************************************************80

void sandia_sgmgg_neighbor_naive ( int dim_num, int point_num, int sparse_index[],
  int neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGG_NEIGHBOR_NAIVE returns the immediate forward neighbor array.
//
//  Discussion:
//
//    A sparse grid index vector is a list of DIM_NUM nonnegative indices.
//
//    An immediate forward L-neighbor of a sparse grid index vector I is a
//    sparse grid index vector J with the property that all entries of J
//    are equal to those of I except for the L-the entry, which is greater by 1.
//
//    A forward neighbor of a sparse grid index vector I is a sparse
//    grid index vector K with the property that every entry of K is
//    equal to or greater by 1 than the corresponding entry of I.
//
//    This function is given a collection of sparse grid index vectors,
//    and returns information defining, for every such vector, the entire
//    set of its immediate forward neighbors.  This is done with a
//    "NEIGHBOR" array of dimension DIM_NUM.  For sparse grid vector I,
//    entry L of NEIGHBOR is 1 if I has an immediate forward L-neighbor,
//    and 0 otherwise.
//
//    This implementation of the procedure is inefficient, and is provided
//    solely for demonstration on small problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
//    Differential Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2411-2442.
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the vector.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int SPARSE_INDEX[DIM_NUM*POINT_NUM],
//    the indices that define the points.
//
//    Output, int NEIGHBOR[DIM_NUM*POINT_NUM], the immediate forward
//    neighbor array.
//
{
  int i;
  int j;
  int j1;
  int j2;
  int l;

  for ( j = 0; j < point_num; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      neighbor[i+j*dim_num] = 0;
    }
  }
  for ( j1 = 0; j1 < point_num; j1++ )
  {
    for ( j2 = 0; j2 < point_num; j2++ )
    {
      l = -1;
      for ( i = 0; i < dim_num; i++ )
      {
//
//  If the entries are not equal...
//
        if ( sparse_index[i+j2*dim_num] != sparse_index[i+j1*dim_num] )
        {
//
//  ...and we haven't already found a difference...
//
          if ( l != -1 )
          {
            l = -1;
            break;
          }
//
//  ...and this difference is +1...
//
          if ( sparse_index[i+j2*dim_num] != sparse_index[i+j1*dim_num] + 1 )
          {
            break;
          }
//
//  ...then remember this index.
//
          l = i;
        }
      }
//
//  If a single unit difference was found, record the direction.
//
      if ( l != - 1 )
      {
        neighbor[l+j1*dim_num] = 1;
      }
    }
  }
  return;
}

}
