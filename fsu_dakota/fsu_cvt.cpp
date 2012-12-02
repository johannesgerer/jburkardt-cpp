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

void fsu_cvt ( int dim_num, int n, int batch, int init, int sample, 
  int sample_num, int it_max, int *seed, double r[], int *it_num )

//****************************************************************************80
//
//  Purpose:
//
//    FSU_CVT computes a Centroidal Voronoi Tessellation.
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
//    This routine initializes the data, and carries out the
//    CVT iteration.
//
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Vance Faber, Max Gunzburger,
//    Centroidal Voronoi Tessellations: Applications and Algorithms,
//    SIAM Review, Volume 41, 1999, pages 637-676.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input, int BATCH, sets the maximum number of sample points
//    generated at one time.  It is inefficient to generate the sample
//    points 1 at a time, but memory intensive to generate them all
//    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
//    BATCH must be at least 1.
//
//    Input, int INIT, specifies how the points are to be initialized.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, points are already initialized on input.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Input/output, int *SEED, the random number seed.
//
//    Input/output, double R[DIM_NUM*N], the approximate CVT points.
//    If INIT = 3 on input, then it is assumed that these values have been
//    initialized.  On output, the CVT iteration has been applied to improve
//    the value of the points.
//
//    Output, int *IT_NUM, the number of iterations taken.  Generally,
//    this will be equal to IT_MAX, unless the iteration tolerance was
//    satisfied early.
//
{
  double energy;
  int i;
  bool reset;
//
//  Check some of the input quantities.
//
  if ( batch < 1 )
  {
    cout << "\n";
    cout << "FSU_CVT - Fatal error!\n";
    cout << "  Input value BATCH < 1.\n";
    exit ( 1 );
  }

  if ( *seed <= 0 )
  {
    cout << "\n";
    cout << "FSU_CVT - Fatal error!\n";
    cout << "  Input value SEED <= 0.\n";
    exit ( 1 );
  }
//
//  Initialize the data.
//
  if ( init != 3 )
  {
    reset = true;

    cvt_sample ( dim_num, n, n, init, reset, seed, r );
  }
//
//  Carry out the iteration.
//
  if ( init == sample )
  {
    reset = false;
  }
  else
  {
    reset = true;
  }

  *it_num = 0;

  for ( i = 1; i <= it_max; i++ )
  {
    *it_num = *it_num + 1;

    cvt_iterate ( dim_num, n, batch, sample, reset, sample_num, seed,
      r, &energy );

    reset = false;
  }
  return;
}
//****************************************************************************80

void cvt_iterate ( int dim_num, int n, int batch, int sample, bool reset, 
  int sample_num, int *seed, double r[], double *energy )

//****************************************************************************80
//
//  Purpose:
//
//    CVT_ITERATE takes one step of the CVT iteration.
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
//    The routine is given a set of points, called "generators", which
//    define a tessellation of the region into Voronoi cells.  Each point
//    defines a cell.  Each cell, in turn, has a centroid, but it is
//    unlikely that the centroid and the generator coincide.
//
//    Each time this CVT iteration is carried out, an attempt is made
//    to modify the generators in such a way that they are closer and
//    closer to being the centroids of the Voronoi cells they generate.
//
//    A large number of sample points are generated, and the nearest generator
//    is determined.  A count is kept of how many points were nearest to each
//    generator.  Once the sampling is completed, the location of all the
//    generators is adjusted.  This step should decrease the discrepancy
//    between the generators and the centroids.
//
//    The centroidal Voronoi tessellation minimizes the "energy",
//    defined to be the integral, over the region, of the square of
//    the distance between each point in the region and its nearest generator.
//    The sampling technique supplies a discrete estimate of this
//    energy.
//
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Vance Faber, and Max Gunzburger,
//    Centroidal Voronoi Tessellations: Applications and Algorithms,
//    SIAM Review, Volume 41, 1999, pages 637-676.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input, int BATCH, sets the maximum number of sample points
//    generated at one time.  It is inefficient to generate the sample
//    points 1 at a time, but memory intensive to generate them all
//    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
//    BATCH must be at least 1.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//
//    Input, bool RESET, is TRUE if the SEED must be reset to SEED_INIT
//    before computation.  Also, the pseudorandom process may need to be
//    reinitialized.
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input/output, int *SEED, the random number seed.
//
//    Input/output, double R[DIM_NUM*N], the Voronoi
//    cell generators.  On output, these have been modified
//
//    Output, double *ENERGY,  the discrete "energy", divided
//    by the number of sample points.
//
{
  int *count;
  int get;
  int have;
  int i;
  int j;
  int j2;
  int *nearest;
  double *r2;
  double *s;
  bool success;
  double term;
//
//  Take each generator as the first sample point for its region.
//  This can slightly slow the convergence, but it simplifies the
//  algorithm by guaranteeing that no region is completely missed
//  by the sampling.
//
  *energy = 0.0;
  r2 = new double[dim_num*n];
  count = new int[n];
  nearest = new int[batch];
  s = new double[dim_num*batch];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      r2[i+j*dim_num] = r[i+j*dim_num];
    }
  }
  for ( j = 0; j < n; j++ )
  {
    count[j] = 1;
  }
//
//  Generate the sampling points S.
//
  have = 0;

  while ( have < sample_num )
  {
    get = i4_min ( sample_num - have, batch );

    cvt_sample ( dim_num, sample_num, get, sample, reset, seed, s );

    reset = false;
    have = have + get;
//
//  Find the index N of the nearest cell generator to each sample point S.
//
    find_closest ( dim_num, n, get, s, r, nearest );
//
//  Add S to the centroid associated with generator N.
//
    for ( j = 0; j < get; j++ )
    {
      j2 = nearest[j];
      for ( i = 0; i < dim_num; i++ )
      {
        r2[i+j2*dim_num] = r2[i+j2*dim_num] + s[i+j*dim_num];
      }
      for ( i = 0; i < dim_num; i++ )
      {
        *energy = *energy + pow ( r[i+j2*dim_num] - s[i+j*dim_num], 2 );
      }
      count[j2] = count[j2] + 1;
    }
  }
//
//  Estimate the centroids.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      r2[i+j*dim_num] = r2[i+j*dim_num] / ( double ) ( count[j] );
    }
  }
//
//  Replace the generators by the centroids.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      r[i+j*dim_num] = r2[i+j*dim_num];
    }
  }
//
//  Normalize the discrete energy estimate.
//
  *energy = *energy / sample_num;

  delete [] count;
  delete [] nearest;
  delete [] r2;
  delete [] s;

  return;
}
//****************************************************************************80

void cvt_sample ( int dim_num, int n, int n_now, int sample, bool reset, 
  int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    CVT_SAMPLE returns sample points.
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
//    N sample points are to be taken from the unit box of dimension DIM_NUM.
//
//    These sample points are usually created by a pseudorandom process
//    for which the points are essentially indexed by a quantity called
//    SEED.  To get N sample points, we generate values with indices
//    SEED through SEED+N-1.
//
//    It may not be practical to generate all the sample points in a 
//    single call.  For that reason, the routine allows the user to
//    request a total of N points, but to require that only N_NOW be
//    generated now (on this call).  
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
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input, int N_NOW, the number of sample points to be generated
//    on this call.  N_NOW must be at least 1.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//
//    Input, bool RESET, is TRUE if the pseudorandom process should be
//    reinitialized.
//
//    Input/output, int *SEED, the random number seed.
//
//    Output, double R[DIM_NUM*N_NOW], the sample points.
//
{
  double exponent;
  static int *halton_base = NULL;
  static int *halton_leap = NULL;
  static int *halton_seed = NULL;
  int halton_step;
  int i;
  int j;
  int k;
  static int ngrid;
  static int rank;
  int rank_max;
  static int *tuple = NULL;

  if ( n_now < 1 )
  {
    cout << "\n";
    cout << "CVT_SAMPLE - Fatal error!\n";
    cout << "  N_NOW < 1.\n";
    exit ( 1 );
  }

  if ( sample == -1 )
  {
    if ( reset )
    {
      random_initialize ( *seed );
    }

    for ( j = 0; j < n_now; j++ )
    {
      for ( i = 0; i < dim_num; i++ )
      {
# ifdef TFLOPS
        r[i+j*dim_num] = ( double ) rand ( ) / ( double ) RAND_MAX;
# else
        r[i+j*dim_num] = ( double ) random ( ) / ( double ) RAND_MAX;
# endif
      }
    }
    *seed = *seed + n_now * dim_num;
  }
  else if ( sample == 0 )
  {
    r8mat_uniform_01 ( dim_num, n_now, seed, r );
  }
  else if ( sample == 1 )
  {
    halton_seed = new int[dim_num];
    halton_leap = new int[dim_num];
    halton_base = new int[dim_num];
      
    halton_step = *seed;

    for ( i = 0; i < dim_num; i++ )
    {
      halton_seed[i] = 0;
    }

    for ( i = 0; i < dim_num; i++ )
    {
      halton_leap[i] = 1;
    }

    for ( i = 0; i < dim_num; i++ )
    {
      halton_base[i] = prime ( i + 1 );
    }

    fsu_halton ( dim_num, n_now, halton_step, halton_seed, halton_leap,
      halton_base, r );

    delete [] halton_seed;
    delete [] halton_leap;
    delete [] halton_base;

    *seed = *seed + n_now;
  }
  else if ( sample == 2 )
  {
    tuple = new int[dim_num];

    exponent = 1.0 / ( double ) ( dim_num );
    ngrid = ( int ) pow ( ( double ) n, exponent );
    rank_max = ( int ) pow ( ( double ) ngrid, ( double ) dim_num );

    if ( rank_max < n )
    {
      ngrid = ngrid + 1;
      rank_max = ( int ) pow ( ( double ) ngrid, ( double ) dim_num );
    }

    if ( reset )
    {
      rank = -1;
      tuple_next_fast ( ngrid, dim_num, rank, tuple );
    }

    rank = ( *seed ) % rank_max;

    for ( j = 0; j < n_now; j++ )
    {
      tuple_next_fast ( ngrid, dim_num, rank, tuple );
      rank = rank + 1;
      rank = rank % rank_max;
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = double ( 2 * tuple[i] - 1 ) / double ( 2 * ngrid );
      }
    }
    *seed = *seed + n_now;
    delete [] tuple;
  }
  else
  {
    cout << "\n";
    cout << "CVT_SAMPLE - Fatal error!\n";
    cout << "  The value of SAMPLE = " << sample << " is illegal.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void cvt_write ( int dim_num, int n, int batch, int seed_init, int seed, 
  char *init_string, int it_max, int it_num, char *sample_string, 
  int sample_num, double r[], char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    CVT_WRITE writes a CVT dataset to a file.
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
//    The initial lines of the file are comments, which begin with a
//    "#" character.
//
//    Thereafter, each line of the file contains the M-dimensional
//    components of the next entry of the dataset.
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
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int BATCH, sets the maximum number of sample points
//    generated at one time.  It is inefficient to generate the sample
//    points 1 at a time, but memory intensive to generate them all
//    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
//
//    Input, int SEED_INIT, the initial random number seed.
//
//    Input, int SEED, the current random number seed.
//
//    Input, char *INIT_STRING, specifies how the initial
//    generators are chosen:
//    filename, by reading data from a file;
//    'GRID', picking points from a grid;
//    'HALTON', from a Halton sequence;
//    'RANDOM', using the C++ RANDOM function;
//    'UNIFORM', using a simple uniform RNG;
//
//    Input, int IT_MAX, the maximum number of iterations allowed.
//
//    Input, int IT_NUM, the actual number of iterations taken.
//
//    Input, char *SAMPLE_STRING, specifies how the region is sampled:
//    'GRID', picking points from a grid;
//    'HALTON', from a Halton sequence;
//    'RANDOM', using the C++ RANDOM function;
//    'UNIFORM', using a simple uniform RNG;
//
//    Input, int SAMPLE_NUM, the number of sampling points used on
//    each iteration.
//
//    Input, double R(DIM_NUM,N), the points.
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
    cout << "CVT_WRITE - Fatal error!\n";
    cout << "  Could not open the output file.\n";
    exit ( 1 );
  }

  s = timestring ( );

  file_out << "#  " << file_out_name << "\n";
  file_out << "#  created by routine CVT_WRITE.CC" << "\n";
  file_out << "#  at " << s << "\n";
  file_out << "#\n";

  file_out << "#  Spatial dimension DIM_NUM =  "  << dim_num          << "\n";
  file_out << "#  Number of points N =      "  << n             << "\n";
  file_out << "#  Initial SEED_INIT =       "  << seed_init     << "\n";
  file_out << "#  Current SEED =            "  << seed          << "\n";
  file_out << "#  INIT =                   \"" << init_string   << "\".\n";
  file_out << "#  Max iterations IT_MAX =   "  << it_max        << "\n";
  file_out << "#  Iterations IT_NUM =       "  << it_num        << "\n";
  file_out << "#  SAMPLE =                 \"" << sample_string << "\".\n";
  file_out << "#  Samples SAMPLE_NUM =      "  << sample_num    << "\n";
  file_out << "#  Sampling BATCH size =     "  << batch         << "\n";
  file_out << "#  EPSILON (unit roundoff) = "  << r8_epsilon ( ) << "\n";

  file_out << "#\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      file_out << setw(10) << r[i+j*dim_num] << "  ";
    }
    file_out << "\n";
  }

  file_out.close ( );

  delete [] s;

  return;
}
