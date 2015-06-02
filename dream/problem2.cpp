# include <cstdlib>
# include <iostream>
# include <cstring>
# include <cmath>

using namespace std;

# include "dream.hpp"
# include "dream_user.hpp"
# include "pdflib.hpp"
# include "problem1_covariance.hpp"

//****************************************************************************80

void problem_size ( int &chain_num, int &cr_num, int &gen_num, int &pair_num, 
  int &par_num )

//****************************************************************************80
/*
  Purpose:

    PROBLEM_SIZE sets information having to do with dimensions.

  Discussion:

    Although this problem involves multivariate normal distributions for
    the initial sample and the density function, in both cases the covariance
    matrix is diagonal, so we don't have to work very hard to set up the
    calculations involving the distributions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2013

  Author:

    John Burkardt

  Reference:

    Jasper Vrugt, CJF ter Braak, CGH Diks, Bruce Robinson, James Hyman, 
    Dave Higdon,
    Accelerating Markov Chain Monte Carlo Simulation by Differential 
    Evolution with Self-Adaptive Randomized Subspace Sampling,
    International Journal of Nonlinear Sciences and Numerical Simulation,
    Volume 10, Number 3, March 2009, pages 271-288.

  Parameters:

    Output, int &CHAIN_NUM, the total number of chains.
    3 <= CHAIN_NUM.

    Output, int &CR_NUM, the total number of CR values.
    1 <= CR_NUM.

    Output, int &GEN_NUM, the total number of generations.
    2 <= GEN_NUM.

    Output, int &PAIR_NUM, the number of pairs of 
    crossover chains.
    0 <= PAIR_NUM.

    Output, int &PAR_NUM, the total number of parameters.
    1 <= PAR_NUM.
*/
{
  chain_num = 10;
  cr_num = 3;
  gen_num = 10;
  pair_num = 3;
  par_num = 10;

  return;
}
//****************************************************************************80

void problem_value ( char **chain_filename, char **gr_filename, 
  double *gr_threshold, int *jumpstep, double limits[], int par_num,
  int *printstep, char **restart_read_filename, char **restart_write_filename )

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
//    01 July 2013
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

  *chain_filename = "problem2_chain00.txt";
  *gr_filename = "problem1_gr.txt";
  gr_threshold = 1.2;
  jumpstep = 5;
  for ( j = 0; j < par_num; j++ )
  {
    limits[0+j*2] = -10.0;
    limits[1+j*2] = +10.0;
  }
  printstep = 10;
  *restart_read_filename = "";
  *restart_write_filename = "problem2_restart.txt";

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
//  Discussion:
//
//    "The initial sample was generated from a normal distribution with
//    variance-covariance matrix 5 * Id."
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2013
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
  double value;
  static double zp_mean_1d = 0.0;
  static double zp_sd_1d = sqrt ( 5.0 );

  value = 1.0;
  for ( i = 0; i < par_num; i++ )
  {
    value = value * r8_normal_pdf ( zp_mean_1d, zp_sd_1d, zp[i] );
  }

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
//  Discussion:
//
//    "The initial sample was generated from a normal distribution with
//    variance-covariance matrix 5 * Id."
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2013
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
  int i;
  double *zp;
  static double zp_mean_1d = 0.0;
  static double zp_sd_1d = sqrt ( 5.0 );

  zp = new double[par_num];

  for ( i = 0; i < par_num; i++ )
  {
    zp(i) = r8_normal_sample ( zp_av_1d, zp_sd_1d );
  }

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
//    "A 10-dimensional twisted Gaussian density function, given by the
//    unnormalized density
//      pi_b(x) proportional to pi(phi_b(x)),
//    with
//      phi_b(x) = (x1, x2+bx1^2-100b,x3,...,x10).
//    Here, pi signifies the density of a multivariate normal distribution,
//    Nd(0,Sigma), with Sigma=diag(100,1,...,1), and phi_b is a function that
//    is used to transform pi to a twisted distribution."
//
//    The value b = 0.01 corresponds to a mildly nonlinear problem,
//    and b = 0.1 to a highly nonlinear problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2013
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
  double b = 0.01;
  int i;
  const double pi = 3.141592653589793;
  double value;
  double *x;
  double xcx;
  double y2;

  y2 = zp[1] + b * zp[0] * zp[0] - 100.0 * b;

  xcx = zp[0] * zp[0] / 100.0 + y2 * y2;
  for ( i = 2; i < par_num; i++ )
  {
    xcx = xcx + zp[i] * zp[i];
  }

  value = - 0.5 * ( double ) ( par_num ) * log ( 2.0 * pi ) 
    - 0.5 * log ( 100.0 ) 
    - 0.5 * xcx;
    
  return value;
}

