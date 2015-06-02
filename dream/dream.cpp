# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "dream.hpp"
# include "dream_user.hpp"
# include "pdflib.hpp"
# include "rnglib.hpp"

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DREAM.
//
//  Discussion:
//
//    The DREAM program was originally developed by Guannan Zhang, of
//    Oak Ridge National Laboratory (ORNL); it has been incorporated into 
//    the DAKOTA package of Sandia National Laboratory, and is
//    intended to form part of the ORNL package known as TASMANIA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2013
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
//  Local parameters:
//
//    Input, string CHAIN_FILENAME, the "base" filename
//    to be used for the chain files.  If this is ""
//    then the chain files will not be written.  This name should 
//    include a string of 0's which will be replaced by the chain 
//    indices.  For example, "chain000.txt" would work as long as the
//    number of chains was 1000 or less.
//
//    Local, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Local, int CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Local, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
//    each sample.
//
//    Local, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Local, double GR[PAR_NUM*GR_NUM], 
//    the Gelman-Rubin R statistic.
//
//    Local, logical GR_CONV, the Gelman-Rubin convergence flag.
//
//    Local, int GR_COUNT, counts the number of generations
//    at which the Gelman-Rubin statistic has been computed.
//
//    Local, string GR_FILENAME, the name of the file
//    in which values of the Gelman-Rubin statistic will be recorded,
//    or the empty string "" if this file is not to be written.
//
//    Local, int GR_NUM, the number of times the Gelman-Rubin
//    statistic may be computed.
//
//    Local, double GR_THRESHOLD, the convergence tolerance for
//    the Gelman-Rubin statistic.
//
//    Local, double JUMPRATE_TABLE[PAR_NUM], the jumprate table.
//
//    Local, int JUMPSTEP, forces a "long jump" every
//    JUMPSTEP generations.
//
//    Local, double LIMITS[2*PAR_NUM], lower and upper bounds
//    for each parameter.
//
//    Local, int PAIR_NUM, the number of pairs of 
//    crossover chains.
//    0 <= PAIR_NUM.
//
//    Local, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Local, int PRINTSTEP, the interval between generations on 
//    which the Gelman-Rubin statistic will be computed and written to a file.
//
//    Local, string RESTART_READ_FILENAME, the name of the file
//    containing restart information.  If the calculation is not a restart,
//    then this should be set to "".
//
//    Local, string RESTART_WRITE_FILENAME, the name of the file
//    to be written, containing restart information.  If a restart file is
//    not to be written, this should be set to "".
//
//    Local, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
{
  string chain_filename = "";
  int chain_num;
  int cr_num;
  double *fit;
  int gen_num;
  double *gr;
  int gr_conv;
  int gr_count;
  string gr_filename = "";
  int gr_num;
  double gr_threshold;
  double *jumprate_table;
  int jumpstep;
  double *limits;
  int pair_num;
  int par_num;
  int printstep;
  string restart_read_filename = "";
  string restart_write_filename = "";
  double *z;

  timestamp ( );

  cout << "\n";
  cout << "DREAM\n";
  cout << "  C++ version\n";
  cout << "  MCMC acceleration by Differential Evolution.\n";
//
//  Get the problem sizes.
//
  problem_size ( chain_num, cr_num, gen_num, pair_num, par_num );
//
//  Decide if the problem sizes are acceptable.
//
  if ( chain_num < 3 )
  {
    cerr << "\n";
    cerr << "DREAM - Fatal error!\n";
    cerr << "  CHAIN_NUM < 1.\n";
    exit ( 1 );
  }
  if ( cr_num < 1 )
  {
    cerr << "\n";
    cerr << "DREAM - Fatal error!\n";
    cerr << "  CR_NUM < 1.\n";
    exit ( 1 );
  }
  if ( gen_num < 2 )
  {
    cerr << "\n";
    cerr << "DREAM - Fatal error!\n";
    cerr << "  GEN_NUM < 2.\n";
    exit ( 1 );
  }
  if ( pair_num < 0 )
  {
    cerr << "\n";
    cerr << "DREAM - Fatal error!\n";
    cerr << "  PAIR_NUM < 0.\n";
    exit ( 1 );
  }
  if ( par_num < 1 )
  {
    cerr << "\n";
    cerr << "DREAM - Fatal error!\n";
    cerr << "  PAR_NUM < 1.\n";
    exit ( 1 );
  }
//
//  Get the problem data values;
//
  limits = r8mat_zero_new ( 2, par_num );

  problem_value ( &chain_filename, &gr_filename, gr_threshold, jumpstep, 
    limits, par_num, printstep, &restart_read_filename, 
    &restart_write_filename );
//
//  Print the data as a job record.
//
  input_print ( chain_filename, chain_num, cr_num, gr_filename, gr_threshold, 
    jumpstep, limits, gen_num, pair_num, par_num, 
    printstep, restart_read_filename, restart_write_filename );
//
//  Allocate and zero out memory.
//
  gr_num = gen_num / printstep;

  fit = r8mat_zero_new ( chain_num, gen_num );
  gr = r8mat_zero_new ( par_num, gr_num );
  z = r8block_zero_new ( par_num, chain_num, gen_num );
//
//  Set the jump rate table.
//
  jumprate_table = jumprate_table_init ( pair_num, par_num );

  jumprate_table_print ( jumprate_table, pair_num, par_num );
//
//  Initialize the Gelman-Rubin data.
//
  gr_init ( gr, gr_conv, gr_count, gr_num, par_num );

  cout << "\n";
  cout << "GR_PRINT:\n";
  cout << "  GR_CONV  = " << gr_conv << "\n";
  cout << "  GR_COUNT = " << gr_count << "\n";
  cout << "  GR_NUM   = " << gr_num << "\n";
//
//  Set the first generation of the chains from restart data, or by sampling.
//
  if ( restart_read_filename.empty ( ) )
  {
    chain_init ( chain_num, fit, gen_num, par_num, z );
  }
  else
  {
    restart_read ( chain_num, fit, gen_num, par_num, restart_read_filename, z );
  }
  chain_init_print ( chain_num, fit, gen_num, par_num, restart_read_filename,
    z );
//
//  Carry out the DREAM algorithm.
//
  dream_algm ( chain_num, cr_num, fit, gen_num, gr, gr_conv, gr_count, 
    gr_num, gr_threshold, jumprate_table, jumpstep, limits, pair_num, 
    par_num, printstep, z );
//
//  Save Gelman-Rubin statistic values to a file.
//
  if ( ! gr_filename.empty ( ) )
  {
    gr_write ( gr, gr_filename, gr_num, par_num, printstep );
  }
//
//  Save parameter values for all chains at last generation.
//
  if ( ! restart_write_filename.empty ( ) )
  {
    restart_write ( chain_num, fit, gen_num, par_num, restart_write_filename,
      z );
  }
//
//  Write each chain to a separate file.
//
  if ( ! chain_filename.empty ( ) )
  {
    chain_write ( chain_filename, chain_num, fit, gen_num, par_num, z );
  }
//
//  Free memory.
//
  delete [] fit;
  delete [] gr;
  delete [] jumprate_table;
  delete [] limits;
  delete [] z;
//
//  Terminate.
//
  cout << "\n";
  cout << "DREAM\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void chain_init ( int chain_num, double fit[], int gen_num, int par_num, 
  double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHAIN_INIT starts Markov chains from a prior distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Output, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
//    each sample.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Output, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
{
  int c;
  int p;
  double *zp;

  for ( c = 0; c < chain_num; c++ )
  {
    zp = prior_sample ( par_num );

    for ( p = 0; p < par_num; p++ )
    {
      z[p+c*par_num+0*par_num*chain_num] = zp[p];
    }

    fit[c+0*chain_num] = sample_likelihood ( par_num, zp );

    delete [] zp;
  }
  return;
}
//****************************************************************************80

void chain_init_print ( int chain_num, double fit[], int gen_num, int par_num, 
  string restart_read_filename, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHAIN_INIT_PRINT prints the initial values for Markov chains.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
//    each sample.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, string RESTART_READ_FILENAME, the name of the 
//    restart file.
//
//    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
{
  int i;
  int j;

  cout << "\n";
  cout << "CHAIN_INIT_PRINT\n";
  cout << "  Display initial values of Markov chains.\n";

  if ( ! restart_read_filename.empty ( ) )
  {
    cout << "  Initialization from restart file \"" 
         << restart_read_filename << "\"\n";
  }
  else
  {
    cout << "  Initialization by sampling prior density.\n";
  }
  for ( j = 0; j < chain_num; j++ )
  {
    cout << "\n";
    cout << "  Chain " << j << "\n";
    cout << "  Fitness " << fit[j+0*chain_num] << "\n";
    for ( i = 0; i < par_num; i++ )
    {
      cout << "  " << setw(14) << z[i+j*par_num+0*par_num*chain_num];
      if ( ( i + 1 ) % 5 == 0 || i == par_num - 1 )
      {
        cout << "\n";
      }
    }
  }
  return;
}
//****************************************************************************80

void chain_outliers ( int chain_num, int gen_index, int gen_num, int par_num,
  double fit[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHAIN_OUTLIERS identifies and modifies outlier chains during burn-in.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
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
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input, int GEN_INDEX, the index of the current generation.
//    2 <= GEN_INDEX <= GEN_NUM.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input/output, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
//    each sample.
//
//    Input/output, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov
//    chain sample data.
//
{
  double *avg;
  double avg_max;
  double *avg_sorted;
  int best;
  int i;
  int ind1;
  int ind3;
  int j;
  int klo;
  int knum;
  int k;
  int outlier_num;
  double q1;
  double q3;
  double qr;
  double t;

  klo = ( ( gen_index + 1 ) / 2 ) - 1;
  knum = gen_index + 1 - klo;

  avg = new double[chain_num];

  for ( j = 0; j < chain_num; j++ )
  {
    t = 0.0;
    for ( k = klo; k <= gen_index; k++ )
    {
      t = t + fit[j+k*chain_num];
    }
    avg[j] = t  / ( double ) ( knum );
  }
//
//  Set BEST to be the index of the chain with maximum average.
//
  best = 0;
  avg_max = avg[0];
  for ( j = 1; j < chain_num; j++ )
  {
    if ( avg_max < avg[j] )
    {
      best = j;
      avg_max = avg[j];
    }
  }
//
//  Determine the indices of the chains having averages 1/4 "above" 
//  and "below" the average.
//
  avg_sorted = r8vec_copy_new ( chain_num, avg );

  r8vec_sort_heap_a ( chain_num, avg_sorted );

  ind1 = r8_round_i4 ( 0.25 * ( double ) ( chain_num ) );
  ind3 = r8_round_i4 ( 0.75 * ( double ) ( chain_num ) );

  q1 = avg_sorted[ind1];
  q3 = avg_sorted[ind3];
  qr = q3 - q1;

  delete [] avg_sorted;
//
//  Identify outlier chains, and replace their later samples
//  with values from the "best" chain.
//
  outlier_num = 0;
  for ( j = 0; j < chain_num; j++ )
  {
    if ( avg[j] < q1 - 2.0 * qr )
    {
      outlier_num = outlier_num + 1;
      for ( i = 0; i < par_num; i++ )
      {
        z[i+j*par_num+gen_index*par_num*chain_num] = 
          z[i+best*par_num+gen_index*par_num*chain_num];
      }
      for ( k = klo; k <= gen_index; k++ )
      {
        fit[j+k*chain_num]  = fit[best+k*chain_num];
      }
    }
  }
//
//  List the outlier chains.
//
  if ( 0 < outlier_num )
  {
    cout << "\n";
    cout << "CHAIN_OUTLIERS:\n";
    cout << "  At iteration " << gen_index 
         << " found " << outlier_num << " outlier chains,\n";
    cout << "  whose indices appear below, and for which samples\n";
    cout << "  from the chain with the largest log likelihood function,\n";
    cout << "  index number " << best << " will be substituted.\n";

    for ( j = 0; j < chain_num; j++ )
    {
      if ( avg[j] < q1 - 2.0 * qr )
      {
        cout << "  " << j << "\n";
      }
    }
  }

  delete [] avg;

  return;
}
//****************************************************************************80

void chain_write ( string chain_filename, int chain_num, double fit[], 
  int gen_num, int par_num, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHAIN_WRITE writes samples of each chain to separate files.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version John Burkardt.
//
//  Parameters:
//
//    Input, string CHAIN_FILENAME, the "base" filename
//    to be used for the chain files.  If this is ""
//    then the chain files will not be written.  This name should 
//    include a string of 0's which will be replaced by the chain 
//    indices.  For example, "chain000.txt" would work as long as the
//    number of chains was 1000 or less.
//
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
//    each sample.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
{
  ofstream chain;
  string chain_filename2;
  int i;
  int j;
  int k;
//
//  Make a temporary copy of the filename template, which we can alter.
//
  chain_filename2 = chain_filename;
//
//  Write parameter samples of all chains.
//
  cout << "\n";
  cout << "CHAIN_WRITE:\n";

  for ( j = 0; j < chain_num; j++ )
  {
    chain.open ( chain_filename2.c_str ( ) );

    if ( ! chain )
    {
      cerr << "\n";
      cerr << "CHAIN_WRITE - Fatal error!\n";
      cerr << "  Could not open file \"" << chain_filename2 << "\".\n";
      exit ( 1 );
    }

    chain <<
      "DREAM.CPP:Parameters_and_log_likelihood_for_chain_#" << j << "\n";

    for ( k = 0; k < gen_num; k++ )
    {
      chain << "  " << k
            << "  " << fit[j+k*chain_num];
      for ( i = 0; i < par_num; i++ )
      {
        chain << "  " << z[i+j*par_num+k*par_num*chain_num];
      }
      chain << "\n";
    }

    chain.close ( );

    cout << "  Created file \"" << chain_filename2 << "\".\n";

    filename_inc ( &chain_filename2 );
  }
  return;
}
//****************************************************************************80

void cr_dis_update ( int chain_index, int chain_num, double cr_dis[], 
  int cr_index, int cr_num, int cr_ups[], int gen_index, int gen_num, 
  int par_num, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    CR_DIS_UPDATE updates the CR distance.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int CHAIN_INDEX, the index of the chain.
//    0 <= CHAIN_INDEX < CHAIN_NUM.
//
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input/output, double CR_DIS[CR_NUM], the CR distances.
//
//    Input, int CR_INDEX, the index of the CR.
//    0 <= CR_INDEX < CR_NUM.
//
//    Input, int CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Input/output, int CR_UPS[CR_NUM], the number of updates 
//    for each CR.
//
//    Input, int GEN_INDEX, the index of the generation.
//    0 <= GEN_INDEX < GEN_NUM.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
{
  int i;
  int i1;
  int i2;
  double *std;
//
//  Compute the standard deviations.
//
  std = std_compute ( chain_num, gen_index, gen_num, par_num, z );
//
//  Increment the update count.
//
  cr_ups[cr_index] = cr_ups[cr_index] + 1;
//
//  Update the CR distance.
//
  for ( i = 0; i < par_num; i++ )
  {
    i1 = i + chain_index * par_num +   gen_index       * par_num * chain_num;
    i2 = i + chain_index * par_num + ( gen_index - 1 ) * par_num * chain_num;
    cr_dis[cr_index] = cr_dis[cr_index] + pow ( ( z[i1] - z[i1] ) / std[i], 2 );
  }

  delete [] std;

  return;
}
//****************************************************************************80

int cr_index_choose ( int cr_num, double cr_prob[] )

//****************************************************************************80
//
//  Purpose:
//
//    CR_INDEX_CHOOSE chooses a CR index.
//
//  Discussion:
//
//    Index I is chosen with probability CR_PROB(I).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Input, double CR_PROB[CR_NUM], the probability of each CR.
//
//    Output, int CR_INDEX_CHOOSE, the index of the CR.
//    0 <= CR_INDEX_CHOOSE < CR_NUM.
//
{
  int cr_index;
  int i;
  int n;
  int *tmp_index;

  if ( cr_num == 1 )
  {
    cr_index = 0;
  }
  else
  { 
    n = 1;
    tmp_index = i4vec_multinomial_sample ( n, cr_prob, cr_num );

    for ( i = 0; i < cr_num; i++ )
    {
      if ( tmp_index[i] == 1 )
      {
        cr_index = i;
        break;
      }
    }
    delete [] tmp_index;
  }
  return cr_index;
}
//****************************************************************************80

void cr_init ( double cr[], double cr_dis[], int cr_num, double cr_prob[], 
  int cr_ups[] )

//****************************************************************************80
//
//  Purpose:
//
//    CR_INIT initializes the crossover probability values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double CR[CR_NUM], the CR values.
//
//    Output, double CR_DIS[CR_NUM], the CR distances.
//
//    Input, int CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Output, double CR_PROB[CR_NUM], the probability of each CR.
//
//    Output, int CR_UPS[CR_NUM], the number of updates
//    for each CR.
//
{
  int i;

  for ( i = 0; i < cr_num; i++ )
  {
    cr[i] = ( double ) ( i + 1 ) / ( double ) ( cr_num );
    cr_dis[i] = 1.0;
    cr_prob[i] = 1.0 / ( double ) ( cr_num );
    cr_ups[i] = 1;
  }
  return;
}
//****************************************************************************80

void cr_prob_update ( double cr_dis[], int cr_num, double cr_prob[], 
  int cr_ups[] )

//****************************************************************************80
//
//  Purpose:
//
//    CR_PROB_UPDATE updates the CR probabilities.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double CR_DIS[CR_NUM], the CR distances.
//
//    Input, int CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Output, double CR_PROB[CR_NUM], the updated CR probabilities.
//
//    Input, int CR_UPS[CR_NUM], the number of updates 
//    for each CR.
//
{
  double cr_prob_sum;
  int i;
 
  for ( i = 0; i < cr_num - 1; i++ )
  {
    cr_prob[i] = cr_dis[i] / ( double ) cr_ups[i];
  }

  cr_prob_sum = r8vec_sum ( cr_num, cr_prob );

  for ( i = 0; i < cr_num - 1; i++ )
  {
    cr_prob[i] = cr_prob[i] / cr_prob_sum;
  }

  return;
}
//****************************************************************************80

double *diff_compute ( int chain_num, int gen_index, int gen_num, 
  int jump_dim[], int jump_num, int pair_num, int par_num, int r[], 
  double z[] ) 

//****************************************************************************80
//
//  Purpose:
//
//    DIFF_COMPUTE computes the differential evolution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
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
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input, int GEN_INDEX, the index of the current generation.
//    1 <= GEN_INDEX <= GEN_NUM.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, int JUMP_DIM[JUMP_NUM], the dimensions in which
//    a jump is to be made.
//
//    Input, int JUMP_NUM, the number of dimensions in which
//    a jump will be made.  0 <= JUMP_NUM <= PAR_NUM.
//
//    Input, int PAIR_NUM, the number of pairs of 
//    crossover chains.
//    0 <= PAIR_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, int R[2*PAIR_NUM], pairs of chains used
//    to compute differences.
//
//    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
//    Output, double DIFF_COMPUTE[PAR_NUM], the vector of pair differences.
//
{
  double *diff;
  int i1;
  int i2;
  int j;
  int k;
  int pair;
  int r1;
  int r2;
//
//  Produce the difference of the pairs used for population evolution.
//
  diff = r8vec_zero_new ( par_num );

  for ( pair = 0; pair < pair_num; pair++ )
  {
    r1 = r[0+pair*2];
    r2 = r[1+pair*2];
    for ( j = 0; j < jump_num; j++ )
    {
      k = jump_dim[j];
      i1 = k+r1*par_num+(gen_index-1)*par_num*chain_num;
      i2 = k+r2*par_num+(gen_index-1)*par_num*chain_num;
      diff[k] = diff[k] + ( z[i1] - z[i2] );
    }
  }

  return diff;
}
//****************************************************************************80

void dream_algm ( int chain_num, int cr_num, double fit[], int gen_num, 
  double gr[], int &gr_conv, int &gr_count, int gr_num, double gr_threshold,
  double jumprate_table[], int jumpstep, double limits[], int pair_num, 
  int par_num, int printstep, double z[] )
    
//****************************************************************************80
//
//  Purpose:
//
//    DREAM_ALGM gets a candidate parameter sample.
//
// Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
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
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input, int CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Input, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
//    each sample.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, double GR[PAR_NUM*GR_NUM], 
//    the Gelman-Rubin R statistic.
//
//    Input/output, int &GR_CONV, the Gelman-Rubin convergence flag.
//
//    Input/output, int &GR_COUNT, counts the number of generations
//    at which the Gelman-Rubin statistic has been computed.
//
//    Input, int GR_NUM, the number of times the Gelman-Rubin
//    statistic may be computed.
//
//    Input, double GR_THRESHOLD, the convergence tolerance for
//    the Gelman-Rubin statistic.
//
//    Input, double JUMPRATE_TABLE[PAR_NUM], the jumprate table.
//
//    Input, int JUMPSTEP, forces a "long jump" every
//    JUMPSTEP generations.
//
//    Input, double LIMITS[2*PAR_NUM], lower and upper bounds
//    for each parameter.
//
//    Input, int PAIR_NUM, the number of pairs of 
//    crossover chains.
//    0 <= PAIR_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, int PRINTSTEP, the interval between generations on 
//    which the Gelman-Rubin statistic will be computed and written to a file.
//
//    Output, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
//  Local parameters:
//
//    Local, int CHAIN_INDEX, the index of the current chain.
//    1 <= CHAIN_INDEX <= CHAIN_NUM.
//
//    Local, double CR[CR_NUM], the CR values.
//
//    Local, double CR_DIS[CR_NUM], the CR distances.
//
//    Local, int CR_INDEX, the index of the selected CR value.
//    1 <= CR_INDEX <= CR_NUM.
//
//    Local, double CR_PROB[CR_NUM], the probability of each CR.
//
//    Local, double CR_UPS[CR_NUM], the number of updates for each CR.
//
//    Local, int GEN_INDEX, the index of the current generation.
//    1 <= GEN_INDEX <= GEN_NUM.
//
//    Local, double ZP[PAR_NUM], a candidate sample.
//
//    Local, int ZP_ACCEPT, the number of candidates accepted.
//
//    Local, double ZP_ACCEPT_RATE, the rate at which generated
//    candidates were accepted.
//
//    Local, int ZP_COUNT, the number of candidates generated.
//
//    Local, double ZP_RATIO, the Metropolis ratio for a candidate.
//
{
  int chain_index;
  double *cr;
  double *cr_dis;
  int cr_index;
  double *cr_prob;
  int *cr_ups;
  int gen_index;
  int i;
  int ind1;
  int ind2;
  double pd1;
  double pd2;
  double r;
  double *zp;
  int zp_accept;
  double zp_accept_rate;
  int zp_count;
  double zp_fit;
  double *zp_old;
  double zp_old_fit;
  double zp_ratio;

  zp_old = new double[par_num];
  zp_count = 0;
  zp_accept = 0;
//
//  Initialize the CR values.
//
  cr = new double[cr_num];
  cr_dis = new double[cr_num];
  cr_prob = new double[cr_num];
  cr_ups = new int[cr_num];

  cr_init ( cr, cr_dis, cr_num, cr_prob, cr_ups );

  for ( gen_index = 1; gen_index < gen_num; gen_index++ )
  {
    for ( chain_index = 0; chain_index < chain_num; chain_index++ )
    {
//
//  Choose CR_INDEX, the index of a CR.
//
      cr_index = cr_index_choose ( cr_num, cr_prob );
//
//  Generate a sample candidate ZP.
//
      zp = sample_candidate ( chain_index, chain_num, cr, cr_index, cr_num, 
        gen_index, gen_num, jumprate_table, jumpstep, limits, pair_num, 
        par_num, z );

      zp_count = zp_count + 1;
//
//  Compute the log likelihood function for ZP.
//
      zp_fit = sample_likelihood ( par_num, zp );

      for ( i = 0; i < par_num; i++ )
      {
        zp_old[i] = z[i+chain_index*par_num+(gen_index-1)*par_num*chain_num];
      }
      zp_old_fit = fit[chain_index+(gen_index-1)*chain_num];
//
//  Compute the Metropolis ratio for ZP.
//
      pd1 = prior_density ( par_num, zp );

      pd2 = prior_density ( par_num, 
        z+0+chain_index*par_num+(gen_index-1)*par_num*chain_num );

      zp_ratio = exp ( 
        ( zp_fit     + log ( pd1 ) ) - 
        ( zp_old_fit + log ( pd2 ) ) );

      zp_ratio = r8_min ( zp_ratio, 1.0 );
//
//  Accept the candidate, or copy the value from the previous generation.
//
      r = r8_uniform_01_sample ( );

      if ( r <= zp_ratio )
      {
        for ( i = 0; i < par_num; i++ )
        {
          z[i+chain_index*par_num+gen_index*par_num*chain_num] = zp[i];
        }
        zp_accept = zp_accept + 1;
        fit[chain_index+gen_index*chain_num] = zp_fit;
      }
      else
      {
        for ( i = 0; i < par_num; i++ )
        {
          z[i+chain_index*par_num+gen_index*par_num*chain_num] = zp_old[i];
        }
        fit[chain_index+gen_index*chain_num] = zp_old_fit; 
      }
//
//  Update the CR distance.
//
      if ( ! gr_conv )
      {
        if ( 1 < cr_num )
        {
          cr_dis_update ( chain_index, chain_num, cr_dis, cr_index, 
            cr_num, cr_ups, gen_index, gen_num, par_num, z );
        }
      }

      delete [] zp;
    }
//
//  Update the multinomial distribution of CR.
//
    if ( ! gr_conv )
    {
      if ( 1 < cr_num )
      {
        if ( ( gen_index + 1 ) % 10 == 0 )
        {
          cr_prob_update ( cr_dis, cr_num, cr_prob, cr_ups );
        }
      }
    }
//
//  Every PRINTSTEP interval,
//  * compute the Gelman Rubin R statistic for this generation,
//    and determine if convergence has occurred.
//
    if ( ( gen_index + 1 ) % printstep == 0 )
    {
      gr_compute ( chain_num, gen_index, gen_num, gr, gr_conv, gr_count, 
        gr_num, gr_threshold, par_num, z );
    }
//
//  Check for outlier chains.
//
    if ( ! gr_conv )
    {
      if ( ( gen_index + 1 ) % 10 == 0 )
      {
        chain_outliers ( chain_num, gen_index, gen_num, par_num, fit, z );
      }
    }
  }    
//
//  Compute the acceptance rate.
//
  zp_accept_rate = ( double ) ( zp_accept ) / ( double ) ( zp_count );

  cout << "\n";
  cout << "  The acceptance rate is " << zp_accept_rate << "\n";

  delete [] cr;
  delete [] cr_dis;
  delete [] cr_prob;
  delete [] cr_ups;
  delete [] zp_old;

  return;
}
//****************************************************************************80

void filename_inc ( string *filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILENAME_INC increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the name is empty, then the routine stops.
//
//    If the name contains no digits, the empty string is returned.
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
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the filename to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILENAME_INC - Fatal error!\n";
    cerr << "  The input string is empty.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens - 1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;

      if ( c == '9' )
      {
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }
//
//  No digits were found.  Return blank.
//
  if ( change == 0 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

  return;
}
//****************************************************************************80

void gr_compute ( int chain_num, int gen_index, int gen_num, double gr[], 
  int &gr_conv, int &gr_count, int gr_num, double gr_threshold, int par_num, 
  double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    GR_COMPUTE computes the Gelman Rubin statistics R used to check
//    convergence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
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
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input, int GEN_INDEX, the index of the current generation.
//    0 < GEN_INDEX < GEN_NUM.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Output, double GR[PAR_NUM*GR_NUM], the Gelman-Rubin R statistic.
//
//    Output, int &GR_CONV, the Gelman-Rubin convergence flag.
//
//    Input/output, int &GR_COUNT, counts the number of 
//    generations at which the Gelman-Rubin statistic has been computed.
//
//    Input, int GR_NUM, the number of times the Gelman-Rubin
//    statistic may be computed.
//
//    Input, double GR_THRESHOLD, the convergence tolerance for the
//    Gelman-Rubin statistic.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
{
  double b_var;
  int chain_index;
  int ind0;
  int k;
  double mean_all;
  double *mean_chain;
  int par_index;
  double rnd0;
  double s;
  double s_sum;
  double var;
  double w_var;

  ind0 = ( ( gen_index + 1 ) / 2 ) - 1;
  rnd0 = ( double ) ( ind0 + 1 );

  mean_chain = new double[chain_num];

  for ( par_index = 0; par_index < par_num; par_index++ )
  {
    for ( chain_index = 0; chain_index < chain_num; chain_index++ )
    {
      mean_chain[chain_index] = 0.0;
      for ( k = ind0; k <= gen_index; k++ )
      {
        mean_chain[chain_index] = mean_chain[chain_index] 
          + z[par_index+chain_index*par_num+k*par_num*chain_num];
      }
      mean_chain[chain_index] = mean_chain[chain_index] / rnd0;
    }

    mean_all = r8vec_sum ( chain_num, mean_chain ) / ( double ) chain_num;

    b_var = 0.0;
    for ( chain_index = 0; chain_index < chain_num; chain_index++ )
    {
      b_var = b_var + pow ( mean_chain[chain_index] - mean_all, 2 );
    }
    b_var = rnd0 * b_var / ( double ) ( chain_num - 1 );

    s_sum = 0.0;
    for ( chain_index = 0; chain_index < chain_num; chain_index++ )
    {
      s = 0.0;
      for ( k = ind0; k <= gen_index; k++ )
      {
        s = s + pow ( z[par_index+chain_index*par_num+k*par_num*chain_num] 
          - mean_chain[chain_index], 2 );
      }
      s_sum = s_sum + s;
    }
    s_sum = s_sum / ( rnd0 - 1.0 );

    w_var = s_sum / ( double ) ( chain_num );

    var = ( ( rnd0 - 1.0 ) * w_var + b_var ) / rnd0;

    gr[par_index+gr_count*par_num] = sqrt ( var / w_var );
  }
//
//  Set the convergence flag.
//
  gr_conv = 1;

  for ( par_index = 0; par_index < par_num; par_index++ )
  {
    if ( gr_threshold < gr[par_index+gr_count*par_num] )
    {
      gr_conv = 0;
      break;
    }
  }

  if ( gr_conv ) 
  {
    cout << "\n";
    cout << "GR_COMPUTE:\n";
    cout << "  GR convergence at iteration: " << gen_index << "\n";
  }

  delete [] mean_chain;

  gr_count = gr_count + 1;

  return;
}
//****************************************************************************80

void gr_init ( double gr[], int &gr_conv, int &gr_count, int gr_num, 
  int par_num )

//****************************************************************************80
//
//  Purpose:
//
//    GR_INIT initializes Gelman-Rubin variables.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double GR[PAR_NUM*GR_NUM], the Gelman-Rubin statistic.
//
//    Output, int &GR_CONV, the convergence flag.
//
//    Output, int &GR_COUNT, counts the number of generations
//    at which the Gelman-Rubin statistic has been computed.
//
//    Input, int GR_NUM, the number of times the Gelman-Rubin
//    statistic may be computed.
//
//    Input, int PAR_NUM, the number of parameters.
//    1 <= PAR_NUM.
//
{
  int i;
  int j;

  for ( j = 0; j < gr_num; j++ )
  {
    for ( i = 0; i < par_num; i++ )
    {
      gr[i+j*par_num] = 0.0;
    }
  }
  gr_conv = 0;
  gr_count = 0;

  return;
}
//****************************************************************************80

void gr_write ( double gr[], string gr_filename, int gr_num, int par_num, 
  int printstep )

//****************************************************************************80
//
//  Purpose:
//
//    GR_WRITE writes Gelman-Rubin R statistics into a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double GR[PAR_NUM*GR_NUM], the Gelman-Rubin R statistic.
//
//    Input, string GR_FILENAME, the Gelman-Rubin filename.
//
//    Input, int GR_NUM, the number of times the Gelman-Rubin
//    statistic may be computed.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, int PRINTSTEP, the interval between generations on 
//    which the Gelman-Rubin statistic will be computed and written to a file.
//
{
  ofstream gr_unit;
  int i;
  int j;

  gr_unit.open ( gr_filename.c_str ( ) );

  if ( ! gr_unit )
  {
    cerr << "\n";
    cerr << "GR_WRITE - Fatal error!\n";
    cerr << "  Could not open the file \"" << gr_filename << "\"\n";
    exit ( 1 );
  }

  gr_unit <<
    "DREAM.CPP:Monitored_parameter_interchains_Gelman_Rubin_R_statistic\n";

  for ( j = 0; j < gr_num; j++ )
  {
    gr_unit << printstep * ( j + 1 ) - 1;
    for ( i = 0; i < par_num; i++ )
    {
      gr_unit << "  " << gr[i+j*par_num];
    }
    gr_unit << "\n";
  }

  gr_unit.close ( );

  cout << "\n";
  cout << "GR_WRITE:\n";
  cout << "  Created the file \"" << gr_filename << "\".\n";

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
//    I4_MIN returns the minimum of two I4's.
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

void i4mat_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT prints an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT_SOME prints some of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << "  " << setw(6) << j - 1;
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to INCX) entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ":";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << "  " << setw(6) << a[i-1+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void i4vec_transpose_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    A = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 }
//    TITLE = "My vector:  "
//
//    My vector:      1    2    3    4    5
//                    6    7    8    9   10
//                   11
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int ihi;
  int ilo;
  int title_len;

  title_len = title.length ( );

  for ( ilo = 1; ilo <= n; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 5 - 1, n );
    if ( ilo == 1 )
    {
      cout << title;
    }
    else
    {
      for ( i = 1; i <= title_len; i++ )
      {
        cout << " ";
      }
    }
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(12) << a[i-1];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

int *i4vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO_NEW creates and zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return a;
}
//****************************************************************************80

void input_print ( string chain_filename, int chain_num, int cr_num, 
  string gr_filename, double gr_threshold, int jumpstep, double limits[], 
  int gen_num, int pair_num, int par_num, int printstep, 
  string restart_read_filename, string restart_write_filename )

//****************************************************************************80
//
//  Purpose:
//
//    INPUT_PRINT prints the data from the input file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CHAIN_FILENAME, the "base" filename
//    to be used for the chain files.  If this is ""
//    then the chain files will not be written.  This name should 
//    include a string of 0's which will be replaced by the chain 
//    indices.  For example, "chain000.txt" would work as long as the
//    number of chains was 1000 or less.
//
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input, int CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Input, string GR_FILENAME, the name of the file
//    in which values of the Gelman-Rubin statistic will be recorded,
//    or the empty string "" if this file is not to be written.
//
//    Input, double GR_THRESHOLD, the convergence tolerance for the
//    Gelman-Rubin statistic.
//
//    Input, int JUMPSTEP, forces a "long jump" every
//    JUMPSTEP generations.
//
//    Input, double LIMITS[2*PAR_NUM], lower and upper limits
//    for each parameter.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, int PAIR_NUM, the number of pairs of 
//    crossover chains.
//    0 <= PAIR_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, int PRINTSTEP, the interval between generations on 
//    which the Gelman-Rubin statistic will be computed and written to a file.
//
//    Local, string RESTART_READ_FILENAME, the name of the file
//    containing restart information.  If the calculation is not a restart,
//    then this should be set to "".
//
//    Local, string RESTART_WRITE_FILENAME, the name of the file
//    to be written, containing restart information.  If a restart file is
//    not to be written, this should be set to "".
//
{
  int j;

  cout << "\n";
  cout << "INPUT_PRINT:\n";
  cout << "\n";
  cout << "  Number of parameters\n";
  cout << "  PAR_NUM = " << par_num << "\n";
  cout << "\n";
  cout << "  LIMITS: Lower and upper limits for each parameter:\n";
  cout << "\n";
  cout << "  Index           Lower           Upper\n";
  cout << "\n";
  for ( j = 0; j < par_num; j++ )
  {
    cout << "  " << setw(5) << j
         << "  " << setw(14) << limits[0+j*2]
         << "  " << setw(14) << limits[1+j*2] << "\n";
  }
  cout << "\n";
  cout << "  Number of generations:\n";
  cout << "  GEN_NUM = " << gen_num << "\n";
  cout << "\n";
  cout << "  Number of simultaneous chains:\n";
  cout << "  CHAIN_NUM = " << chain_num << "\n";
  cout << "\n";
  cout << "  Chain filename (base):\n";
  if ( chain_filename.empty ( ) )
  {
    cout << "  CHAIN_FILENAME = \"(Null)\".\n";
  }
  else
  {
    cout << "  CHAIN_FILENAME = \"" << chain_filename << "\".\n";
  }
  cout << "\n";
  cout << "  Number of pairs of chains for crossover:\n";
  cout << "  PAIR_NUM = " << pair_num << "\n";
  cout << "\n";
  cout << "  Number of crossover values:\n";
  cout << "  CR_NUM = " << cr_num << "\n";
  cout << "\n";
  cout << "  Number of steps til a long jump:\n";
  cout << "  JUMPSTEP = " << jumpstep << "\n";
  cout << "\n";
  cout << "  Interval between Gelman-Rubin computations:\n";
  cout << "  PRINTSTEP = " << printstep << "\n";
  cout << "\n";
  cout << "  Gelman-Rubin data filename:\n";
  if ( gr_filename.empty ( ) )
  {
    cout << "  GR_FILENAME = \"(Null)\".\n";
  }
  else
  {
    cout << "  GR_FILENAME = \"" << gr_filename << "\".\n";
  }
  cout << "  GR_THRESHOLD = " << gr_threshold << "\n";
  cout << "\n";
  cout << "  Restart read filename:\n";
  if ( restart_read_filename.empty ( ) )
  {
    cout << "  RESTART_READ_FILENAME = \"(Null)\".\n";
  }
  else
  {
    cout << "  RESTART_READ_FILENAME = \"" << restart_read_filename << "\".\n";
  }
  cout << "\n";
  cout << "  Restart write filename:\n";
  if ( restart_write_filename.empty ( ) )
  {
    cout << "  RESTART_WRITE_FILENAME = \"(Null)\".\n";
  }
  else
  {
    cout << "  RESTART_WRITE_FILENAME = \"" << restart_write_filename << "\".\n";
  }

  return;
}
//****************************************************************************80

void jumprate_choose ( double cr[], int cr_index, int cr_num, int gen_index,
  int jump_dim[], int &jump_num, double &jumprate, double jumprate_table[],
  int jumpstep, int par_num )

//****************************************************************************80
//
//  Purpose:
//
//    JUMPRATE_CHOOSE chooses a jump rate from the jump rate table.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
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
//    Input, double CR[CR_NUM], the CR values.
//
//    Input, int CR_INDEX, the index of the CR.
//    1 <= CR_INDEX <= CR_NUM.
//
//    Input, int CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Input, int GEN_INDEX, the current generation.
//    1 <= GEN_INDEX <= GEN_NUM.
//
//    Output, int JUMP_DIM[PAR_NUM], the indexes of the
//    parameters to be updated.
//
//    Output, int &JUMP_NUM, the number of dimensions in which
//    a jump will be made.  0 <= JUMP_NUM <= PAR_NUM.
//
//    Output, double &JUMPRATE, the jump rate.
//
//    Input, double JUMPRATE_TABLE[PAR_NUM], the jump rate table.
//
//    Input, int JUMPSTEP, forces a "long jump" every
//    JUMPSTEP generations.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
{
  int i;
  double r;
//
//  Determine the dimensions that will be updated.
//
  jump_num = 0;
  for ( i = 0; i < par_num; i++ )
  {
    jump_dim[i] = 0;
  }

  for ( i = 0; i < par_num; i++ )
  {
    r = r8_uniform_01_sample ( );

    if ( 1.0 - cr[cr_index] < r )
    {
      jump_dim[jump_num] = i;
      jump_num = jump_num + 1;
    }
  }
//
//  Calculate the general jump rate.
//
  if ( jump_num == 0 )
  {
    jumprate = 0.0;
  }
  else
  {
    jumprate = jumprate_table[jump_num-1];
  }
//
//  If parameter dimension is 1, 2, or 3, fix the jump rate to 0.6.
//
  if ( par_num <= 3 )
  {
    jumprate = 0.6;
  }
//
//  Determine if a long jump is forced.
//
  if ( ( gen_index % jumpstep ) == 0 )
  {
    jumprate = 0.98;
  }
  return;
}
//****************************************************************************80

double *jumprate_table_init ( int pair_num, int par_num )

//****************************************************************************80
//
//  Purpose:
//
//    JUMPRATE_TABLE_INIT initializes the jump rate table.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int PAIR_NUM, the number of pairs of 
//    crossover chains.
//    0 <= PAIR_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Output, double JUMPRATE_TABLE_INIT[PAR_NUM], the jumprate table.
//
{
  double c;
  int i;
  double *jumprate_table;
 
  jumprate_table = new double[par_num];

  c = 2.38 / sqrt ( ( double ) ( 2 * pair_num ) );

  for ( i = 0; i < par_num; i++ )
  {
    jumprate_table[i] = c / sqrt ( ( double ) ( i + 1 ) );
  }

  return jumprate_table;
}
//****************************************************************************80

void jumprate_table_print ( double jumprate_table[], int pair_num, int par_num )

//****************************************************************************80
//
//  Purpose:
//
//    JUMPRATE_TABLE_PRINT prints the jump rate table.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double JUMPRATE_TABLE[PAR_NUM], the jumprate table.
//
//    Input, int PAIR_NUM, the number of pairs of 
//    crossover chains.
//    0 <= PAIR_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
{
  int i;
 
  cout << "\n";
  cout << "JUMPRATE_TABLE_PRINT\n";
  cout << "\n";
  cout << "   I        Jumprate\n";
  cout << "\n";
  for ( i = 0; i < par_num; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(14) << jumprate_table[i] << "\n";
  }
  return;
}
//****************************************************************************80

int r8_round_i4 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ROUND_I4 rounds an R8, returning an I4.
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
//    25 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_ROUND_I4, the rounded value.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - floor ( - x + 0.5 );
  }
  else
  {
    value =   floor (   x + 0.5 );
  }

  return value;
}
//****************************************************************************80

double *r8block_zero_new ( int l, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLOCK_ZERO_NEW returns a new zeroed R8BLOCK.
//
//  Discussion:
//
//    An R8BLOCK is a triple dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int L, M, N, the number of rows and columns.
//
//    Output, double R8BLOCK_ZERO_NEW[L*M*N], the new zeroed matrix.
//
{
  double *a;
  int i;
  int j;
  int k;

  a = new double[l*m*n];

  for ( k = 0; k < n; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( i = 0; i < l; i++ )
      {
        a[i+j*l+k*l*m] = 0.0;
      }
    }
  }
  return a;
}
//****************************************************************************80

double *r8mat_zero_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_ZERO_NEW returns a new zeroed R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, double R8MAT_ZERO_NEW[M*N], the new zeroed matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}
//****************************************************************************80

double *r8vec_copy_new ( int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY_NEW copies an R8VEC to a new R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double R8VEC_COPY_NEW[N], the copy of A1.
//
{
  double *a2;
  int i;

  a2 = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
//****************************************************************************80

void r8vec_heap_d ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_HEAP_D reorders an R8VEC into a descending heap.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    A heap is an array A with the property that, for every index J,
//    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
//    2*J+1 and 2*J+2 are legal).
//
//  Diagram:
//
//                  A(0)
//
//            A(1)         A(2)
//
//      A(3)       A(4)  A(5) A(6)
//
//    A(7) A(8)  A(9) A(10)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the input array.
//
//    Input/output, double A[N].
//    On input, an unsorted array.
//    On output, the array has been reordered into a heap.
//
{
  int i;
  int ifree;
  double key;
  int m;
//
//  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
//
  for ( i = (n/2)-1; 0 <= i; i-- )
  {
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = a[i];
    ifree = i;

    for ( ; ; )
    {
//
//  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
//  IFREE.  (One or both may not exist because they equal or exceed N.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n <= m )
      {
        break;
      }
      else
      {
//
//  Does the second position exist?
//
        if ( m + 1 < n )
        {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }
      }
    }
//
//  When you have stopped shifting items up, return the item you
//  pulled out back to the heap.
//
    a[ifree] = key;
  }

  return;
}
//****************************************************************************80

void r8vec_sort_heap_a ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_A ascending sorts an R8VEC using heap sort.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, double A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  int n1;
  double temp;

  if ( n <= 1 )
  {
    return;
  }
//
//  1: Put A into descending heap form.
//
  r8vec_heap_d ( n, a );
//
//  2: Sort A.
//
//  The largest object in the heap is in A[0].
//  Move it to position A[N-1].
//
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
//
//  Consider the diminished heap of size N1.
//
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
//
//  Restore the heap structure of the initial N1 entries of A.
//
    r8vec_heap_d ( n1, a );
//
//  Take the largest object from A[0] and move it to A[N1-1].
//
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }

  return;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}
//****************************************************************************80

void r8vec_transpose_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Example:
//
//    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
//    TITLE = 'My vector:  '
//
//    My vector:
//
//        1.0    2.1    3.2    4.3    5.4
//        6.5    7.6    8.7    9.8   10.9
//       11.0
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
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int ihi;
  int ilo;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= 0 )
  {
    cout << "  (Empty)\n";
    return;
  }

  for ( ilo = 0; ilo < n; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 5, n );
    for ( i = ilo; i < ihi; i++ )
    {
      cout << "  " << setw(12) << a[i];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

double *r8vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO_NEW creates and zeroes an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double R8VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return a;
}
//****************************************************************************80

void restart_read ( int chain_num, double fit[], int gen_num, int par_num, 
  string restart_read_filename, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    RESTART_READ reads parameter sample data from a restart file.
//
//  Discussion:
//
//    Only a single generation (presumably the last one) was written to the file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Output, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
//    each sample.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, string RESTART_READ_FILENAME, the name of 
//    the restart file.
//
//    Output, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
{
  int chain_index;
  int gen_index = 0;
  int index;
  string line;
  size_t n;
  int par_index;
  ifstream restart;

  restart.open ( restart_read_filename.c_str ( ) );

  if ( ! restart )
  {
    cerr << "\n";
    cerr << "RESTART_READ - Fatal error!\n";
    cerr << "  Could not open the file \"" 
         << restart_read_filename << "\".\n";
    exit ( 1 );
  }
//
//  Read and ignore line 1.
//
  getline ( restart, line );
//
//  Read the final fitness and parameter values for each chain.
//
  for ( chain_index = 0; chain_index < chain_num; chain_index++ )
  {
    index = chain_index 
          + chain_num * gen_index;
    restart >> fit[index];
    for ( par_index = 0; par_index < par_num; par_index++ )
    {
      index = par_index 
            + par_num * chain_index 
            + par_num * chain_num * gen_index;
      restart >> z[index];
    }
  }

  restart.close ( );

  return;
}
//****************************************************************************80

void restart_write ( int chain_num, double fit[], int gen_num, int par_num, 
  string restart_write_filename, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    RESTART_WRITE writes a restart file.
//
//  Discussion:
//
//    Only data for the final generation is written.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input, double FIT[CHAIN_NUM*GEN_NUM], the likelihood of
//    each sample.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, string RESTART_WRITE_FILENAME, the name of the 
//    restart file.
//
//    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
{
  int c;
  int p;
  ofstream restart;

  restart.open ( restart_write_filename.c_str ( ) );   

  if ( !restart )
  {
    cerr << "\n";
    cerr << "RESTART_WRITE - Fatal error!\n";
    cerr << "  Could not open \"" << restart_write_filename << "\".\n";
    exit ( 1 );
  }

  restart << "DREAM.C:Parameter_values_for_restart.\n";

  for ( c = 0; c < chain_num; c++ )
  {
    restart << "  " << c
            << "  "<< fit[c+(gen_num-1)*chain_num];
    for ( p = 0; p < par_num; p++ )
    {
      restart << "  " << z[p+c*par_num+(gen_num-1)*par_num*chain_num];
    }
    restart << "\n"; 
  }

  restart.close ( );

  cout << "\n";
  cout << "RESTART_WRITE:\n";
  cout << "  Created restart file \"" << restart_write_filename << "\".\n";

  return;
}
//****************************************************************************80

double *sample_candidate ( int chain_index, int chain_num, double cr[], 
  int cr_index, int cr_num, int gen_index, int gen_num, 
  double jumprate_table[], int jumpstep, double limits[], int pair_num, 
  int par_num, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_CANDIDATE generates candidate parameter samples.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
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
//    Input, int CHAIN_INDEX, the chain index.
//    0 <= CHAIN_INDEX < CHAIN_NUM.
//
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input, double CR[CR_NUM], the CR values.
//
//    Input, int CR_INDEX, the index of the chosen CR value.
//    0 <= CR_INDEX < CR_NUM.
//
//    Input, int CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Input, int GEN_INDEX, the current generation.
//    0 <= GEN_INDEX < GEN_NUM.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, double JUMPRATE_TABLE[PAR_NUM], the jumprate table.
//
//    Input, int JUMPSTEP, forces a "long jump" every
//    JUMPSTEP generations.
//
//    Input, double LIMITS[2*PAR_NUM], limits for the parameters.
//
//    Input, int PAIR_NUM, the number of pairs of 
//    crossover chains.
//    0 <= PAIR_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
//    Output, double SAMPLE_CANDIDATE[PAR_NUM], a candidate parameter sample.
//
//  Local parameters:
//
//    Input, int JUMP_DIM[JUMP_NUM], the dimensions in which
//    a jump is to be made.
//
//    Local, int JUMP_NUM, the number of dimensions in which
//    a jump will be made.  0 <= JUMP_NUM <= PAR_NUM.
//
//    Local, double JUMPRATE, the jump rate.
//
{
  double av;
  double b;
  double *diff;
  double *eps;
  int i;
  int *jump_dim;
  int jump_num;
  double jumprate;
  double *noise_e;
  int pair[2];
  int *r;
  double r2;
  double sd;
  double *zp;
//
//  Used to calculate E following a uniform distribution on (-B,+B).
//  Because B is currently zero, the noise term is suppressed.
//
  b = 0.0;
//
//  Pick pairs of other chains for crossover.
//
  r = new int[2*pair_num];

  for ( i = 0; i < pair_num; i++ )
  {
    while ( 1 )
    {
      r2 = r8_uniform_01_sample ( );
      pair[0] = ( int ) ( r2 * ( double ) chain_num );
      r2 = r8_uniform_01_sample ( );
      pair[1] = ( int ) ( r2 * ( double ) chain_num );

      if ( pair[0] != pair[1] &&
           pair[0] != chain_index && 
           pair[1] != chain_index )
      {
        break;
      }
    }
    r[0+i*2] = pair[0];
    r[1+i*2] = pair[1];
  }
//
//  Determine the jump rate.
//
  jump_dim = new int[par_num];

  jumprate_choose ( cr, cr_index, cr_num, gen_index, jump_dim, jump_num, 
    jumprate, jumprate_table, jumpstep, par_num );
//
//  Calculate E in equation 4 of Vrugt.
//
  noise_e = new double[par_num];

  for ( i = 0; i < par_num; i++ )
  {
    noise_e[i] = b * ( 2.0 * r8_uniform_01_sample ( ) - 1.0 );
  }
//
//  Get epsilon value from multinormal distribution                      
//
  eps = new double[par_num];

  av = 0.0;
  sd = 1.0E-10;
  for ( i = 0; i < par_num; i++ )
  {
    eps[i] = r8_normal_sample ( av, sd );
  }
//
//  Generate the candidate sample ZP based on equation 4 of Vrugt.
//
  diff = diff_compute ( chain_num, gen_index, gen_num, jump_dim, jump_num, 
    pair_num, par_num, r, z );

  zp = new double[par_num];

  for ( i = 0; i < par_num; i++ )
  {
    zp[i] = z[i+chain_index*par_num+(gen_index-1)*par_num*chain_num];
  }
  for ( i = 0; i < par_num; i++ )
  {
    zp[i] = zp[i] + ( 1.0 + noise_e[i] ) * jumprate * diff[i] + eps[i];
  }
//
//  Enforce limits on the sample ZP.
//
  sample_limits ( limits, par_num, zp );

  delete [] diff;
  delete [] eps;
  delete [] jump_dim;
  delete [] noise_e;
  delete [] r;

  return zp;
}
//****************************************************************************80

void sample_limits ( double limits[], int par_num, double zp[] )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_LIMITS enforces limits on a sample variable.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double LIMITS[2*PAR_NUM], the parameter limits.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input/output, double ZP[PAR_NUM], a variable, whose entries,
//    if necessary, will be "folded" so that they lie within the limits.
//
{
  int i;
  double w;

  for ( i = 0; i < par_num; i++ )
  {
    w = limits[1+i*2] - limits[0+i*2];

    if ( w == 0.0 )
    {
      zp[i] = limits[0+i*2];
    }
    else if ( w < 0.0 )
    {
      cerr << "\n";
      cerr << "SAMPLE_LIMITS - Fatal error!\n";
      cerr << "  Upper limit less than lower limit.\n";
      exit ( 1 );
    }
    else
    {
      while ( zp[i] < limits[0+i*2] )
      {
        zp[i] = zp[i] + w;
      }
      while ( limits[1+i*2] < zp[i] )
      {
        zp[i] = zp[i] - w;
      }
    }
  }
  return;
}
//****************************************************************************80

double *std_compute ( int chain_num, int gen_index, int gen_num, int par_num, 
  double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    STD_COMPUTE computes the current standard deviations, for each parameter.
//
//  Discussion:
//
//    The computation encompasses all chains and generations up to the
//    current ones.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Input, int GEN_INDEX, the current generation.
//    0 <= GEN_INDEX < GEN_NUM.
//
//    Input, int GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, double Z[PAR_NUM*CHAIN_NUM*GEN_NUM], the Markov chain 
//    sample data.
//
//    Output, double STD_COMPUTE[PAR_NUM], the standard deviations.
//
{
  int i;
  int j;
  int k;
  double mean;
  double *std;

  std = new double[par_num];

  for ( i = 0; i < par_num; i++ )
  {
    mean = 0.0;
    for ( k = 0; k <= gen_index; k++ )
    {
      for ( j = 0; j < chain_num; j++ )
      {
        mean = mean + z[i+j*par_num+k*par_num*chain_num];
      }
    }
    mean = mean / ( double ) ( chain_num ) / ( double ) ( gen_index );

    std[i] = 0.0;
    for ( k = 0; k <= gen_index; k++ )
    {
      for ( j = 0; j < chain_num; j++ )
      {
        std[i] = std[i] + pow ( z[i+j*par_num+k*par_num*chain_num] - mean, 2 );
      }
    }
    std[i] = sqrt ( std[i] / ( double ) ( chain_num * gen_index - 1 ) );
  }

  return std;
}

