# include <cstdlib>
# include <iostream>
# include <cstring>
# include <cmath>

using namespace std;

# include "dream.hpp"
# include "dream_user.hpp"
# include "pdflib.hpp"
# include "problem1c_covariance.hpp"

Covariance c ( 5 );

//****************************************************************************80

void problem_size ( int &chain_num, int &cr_num, int &gen_num, int &pair_num, 
  int &par_num )

//****************************************************************************80
//
//  Purpose:
//
//    PROBLEM_SIZE sets information having to do with dimensions.
//
//  Discussion:
//
//    In the Vrugt paper, PAR_NUM is 100.  For testing, it is reasonable
//    to try a tiny value like PAR_NUM = 5 instead.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jasper Vrugt, CJF ter Braak, CGH Diks, Bruce Robinson, James Hyman, 
//    Dave Higdon,
//    Accelerating Markov Chain Monte Carlo Simulation by Differential 
//    Evolution with Self-Adaptive Randomized Subspace Sampling,
//    International Journal of Nonlinear Sciences and Numerical Simulation,
//    Volume 10, Number 3, March 2009, pages 271-288.
//
//  Parameters:
//
//    Output, int &CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Output, int &CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Output, int &GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Output, int &PAIR_NUM, the number of pairs of 
//    crossover chains.
//    0 <= PAIR_NUM.
//
//    Output, int &PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
{
  chain_num = 10;
  cr_num = 3;
  gen_num = 10;
  pair_num = 3;
//par_num = 100;
  par_num = 5;

  return;
}
//****************************************************************************80

void problem_value ( string *chain_filename, string *gr_filename, 
  double &gr_threshold, int &jumpstep, double limits[], int par_num, 
  int &printstep, string *restart_read_filename, 
  string *restart_write_filename )

//****************************************************************************80
//
//  Purpose:
//
//    PROBLEM_VALUE sets information, including numeric data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string CHAIN_FILENAME, the "base" filename
//    to be used for the chain files.  If this is ""
//    then the chain files will not be written.  This name should 
//    include a string of 0's which will be replaced by the chain 
//    indices.  For example, "chain000.txt" would work as long as the
//    number of chains was 1000 or less.
//
//    Output, string *GR_FILENAME, the name of the file
//    in which values of the Gelman-Rubin statistic will be recorded,
//    or "" if this file is not to be written.
//
//    Output, double &GR_THRESHOLD, the convergence tolerance for
//    the Gelman-Rubin statistic.
//
//    Output, int &JUMPSTEP, forces a "long jump" every
//    JUMPSTEP generations.
//
//    Output, double LIMITS[2*PAR_NUM], lower and upper bounds
//    for each parameter.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Output, int &PRINTSTEP, the interval between generations on 
//    which the Gelman-Rubin statistic will be computed and written to a file.
//
//    Output, string *RESTART_READ_FILENAME, the name of the file
//    containing restart information.  If this calculation is not a restart,
//    then this should be "".
//
//    Output, string *RESTART_WRITE_FILENAME, the name of the file
//    to be written, containing restart information.  If a restart file is not
//    to be written, this should be "".
//
{
  int j;

  *chain_filename = "problem1_chain00.txt";
  *gr_filename = "problem1_gr.txt";
  gr_threshold = 1.2;
  jumpstep = 5;
  for ( j = 0; j < par_num; j++ )
  {
    limits[0+j*2] =   9.9;
    limits[1+j*2] = +10.0;
  }
  printstep = 10;
  *restart_read_filename = "";
  *restart_write_filename = "problem1_restart.txt";

  return;
}
//****************************************************************************80

double prior_density ( int par_num, double zp[] )

//****************************************************************************80
//
//  Purpose:
//
//    PRIOR_DENSITY evaluates the prior density function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, double ZP[PAR_NUM], the argument of the density
//    function.
//
//    Output, real PRIOR_DENSITY, the value of the prior density function.
//
{
  double det;
  double *factor;
  double *mean;
  double value;

  mean = c.mean_get ( );
  factor = c.factor_get ( );
  det = c.det_get ( );

  value = r8vec_multinormal_pdf ( par_num, mean, factor, det, zp );

  delete [] factor;
  delete [] mean;

  return value;
}
//****************************************************************************80

double *prior_sample ( int par_num )

//****************************************************************************80
//
//  Purpose:
//
//    PRIOR_SAMPLE samples from the prior distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Output, double PRIOR_SAMPLE[PAR_NUM], the sample from the distribution.
//
{
  double *factor;
  int i;
  double *mean;
  double *x;
  double *zp;

  x = new double[par_num];

  for ( i = 0; i < par_num; i++ )
  {
    x[i] = r8_normal_01_sample ( );
  }

  factor = c.factor_get ( );
  zp = r8mat_mtv_new ( par_num, par_num, factor, x );

  mean = c.mean_get ( );
  for ( i = 0; i < par_num; i++ )
  {
    zp[i] = zp[i] + mean[i];
  }

  delete [] factor;
  delete [] mean;
  delete [] x;

  return zp;
}
//****************************************************************************80

double sample_likelihood ( int par_num, double zp[] )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_LIKELIHOOD computes the log likelihood function.
//
//  Discussion:
//
//    This is a one mode Gaussian.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, double ZP[PAR_NUM], a sample.
//
//    Output, double SAMPLE_LIKELIHOOD, the log likelihood function 
//    for the sample.
//
{
  double det;
  double *factor;
  int i;
  double *mean;
  const double pi = 3.141592653589793;
  double value;
  double *x;
  double xcx;
  double *y;

  mean = c.mean_get ( );
  x = new double[par_num];

  for ( i = 0; i < par_num; i++ )
  {
    x[i] = zp[i] - mean[i];
  }

  factor = c.factor_get ( );
  y = r8mat_utsol ( par_num, factor, x );
//
//  Compute:
//    (x-mu)' * inv(C)          * (x-mu)
//  = (x-mu)' * inv(R'*R)       * (x-mu)
//  = (x-mu)' * inv(R) * inv(R) * (x-mu)
//  = y' * y.
//
  xcx = r8vec_dot_product ( par_num, y, y );

  det = c.det_get ( );
  value = - 0.5 * ( double ) ( par_num ) * log ( 2.0 * pi ) 
    - 0.5 * log ( det ) 
    - 0.5 * xcx;

  delete [] factor;
  delete [] mean;
  delete [] x;
  delete [] y;
    
  return value;
}

