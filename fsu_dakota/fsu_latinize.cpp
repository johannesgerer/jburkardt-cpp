# ifdef ANSI_HEADERS
#   include <cstdlib>
#   include <cmath>
#   include <ctime>
# else
#   include <stdlib.h>
#   include <math.h>
#   include <time.h>
# endif

# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "fsu.H"

//****************************************************************************80

void fsu_latinize ( int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    FSU_LATINIZE "Latinizes" a real table dataset.
//
//  License:
//
//    Copyright (C) 2004  John Burkardt and Max Gunzburger
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Discussion:
//
//    It is assumed, though not necessary, that the input dataset
//    has points that lie in the unit hypercube.
//
//    In any case, the output dataset will have this property.
//
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of cells.
//
//    Input/output, double TABLE[M*N].  On input, the dataset to
//    be "Latinized".  On output, the Latinized dataset.
//
{
  double *column;
  int i;
  int *indx;
  int j;

  column = new double[n];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      column[j] = table[i+j*m];
    }
    indx = r8vec_sort_heap_index_a ( n, column );

    for ( j = 0; j < n; j++ )
    {
      table[i+indx[j]*m] = ( double ) ( 2 * j + 1 ) / ( double ) ( 2 * n );
    }

    delete [] indx;
  }

  delete [] column;

  return;
}
