# include "sandia_rules.hpp"
# include "sparse_grid_mixed.hpp"

# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
int rule_string_to_index ( string rule_string );
bool s_eqi ( string s1, string s2 );
void sparse_grid_mixed_dataset_handle ( int dim_num, int level_max, int rule[], 
  double alpha[], double beta[], double tol, std::string file_name );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_MIXED_DATASET.
//
//  Discussion:
//
//    This program computes a sparse grid quadrature rule based on a mixture
//    of 1D rules, and writes it to a file.
//
//    The user specifies:
//    * M, the spatial dimension of the quadrature region,
//    * L, the level that defines the Smolyak grid.
//
//    Then the user specifies rules for each of the M dimensions.
//    A rule, when specified, may be used for one, or for multiple consecutive 
//    dimensions.
//
//    * RULE identifies the 1D rule.
//      "CC", "F2", "GP", "GL", "GH", "GGH", "LG", "GLG", "GJ", "GW",
//      "CCS", "F2S", "GPS".
//    * the number of times the rule is to be used.
//    * ALPHA parameter for that rule;
//    * BETA parameter for that rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2009
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
{
  double *alpha;
  double alpha_1d;
  double *beta;
  double beta_1d;
  int dim;
  int dim_inc;
  int dim_index;
  int dim_num;
  string file_name;
  int level_max;
  int *rule;
  int rule_1d;
  string rule_string;
  double tol;

  webbur::timestamp ( );
  cout << "\n";
  cout << "SPARSE_GRID_MIXED_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Compute the abscissas and weights of a quadrature rule\n";
  cout << "  associated with a sparse grid derived from a Smolyak\n";
  cout << "  construction based on a mixture of 1D rules.\n";
  cout << "\n";
  cout << "  Inputs to the program include:\n";
  cout << "\n";
  cout << "    DIM_NUM, the spatial dimension.\n";
  cout << "    (typically in the range of 2 to 10)\n";
  cout << "\n";
  cout << "    LEVEL_MAX, the level of the sparse grid.\n";
  cout << "    (typically in the range of 0, 1, 2, 3, ...\n";
  cout << "\n";
  cout << "   Then the user must define 1D quadrature rules to be used.\n";
  cout << "   Each rule is used for at least the \"next\" dimension, but can be\n";
  cout << "   used for several or all the remaining consecutive dimensions.\n";
  cout << "\n";
  cout << "  Rule definition requires:\n";
  cout << "\n";
  cout << "  * Rule identifier:\n";
  cout << "    CC, F2, GP, GL, GH, GGH, LG, GLG, GJ, GW, CCS, F2S, GPS.\n";
  cout << "  * Repetition factor (consecutive dimensions with same rule):\n";
  cout << "  * ALPHA, (only for GGH, GLG, GJ rules)\n";
  cout << "  * BETA, (only for GJ rule.)\n";
  cout << "  Output from the program includes:\n";
  cout << "\n";
  cout << "    * Files that define the alphas, betas, ranges, weights and abscissas.\n";
//
//  Get the spatial dimension.
//
  if ( 1 < argc )
  {
    dim_num = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of DIM_NUM (1 or greater)\n";
    cin >> dim_num;
  }

  cout << "\n";
  cout << "  Spatial dimension requested is = " << dim_num << "\n";
//
//  Get the level.
//
  if ( 2 < argc )
  {
    level_max = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of LEVEL_MAX (0 or greater).\n";
    cin >> level_max;
  }
  cout << "  LEVEL_MAX is = " << level_max << "\n";
//
//  Now get the rules.
//
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];

  dim_index = 0;

  while ( dim_index < dim_num )
  {
    cout << "\n";
    cout << "  Rule identifiers include:\n";
    cout << "  CC, F2, GP, GL, GH, GGH, LG, GLG, GJ, GW, CCS, F2S, GPS\n";
    cout << "  Enter the rule identifier for dimension " << dim_index + 1 << "\n";
    cin >> rule_string;
    rule_1d = rule_string_to_index ( rule_string );

    cout << "\n";
    cout << "  How many consecutive dimensions will this same rule be used?\n";
    cin >> dim_inc;

    if ( dim_num < dim_index + dim_inc )
    {
      cout << "\n";
      cout << "SPARSE_GRID_MIXED_DATASET - Fatal error!\n";
      cout << "  Dimension count exceeds limit.\n";
      exit ( 1 );
    }

    if ( rule_1d == 6 || rule_1d == 8 || rule_1d == 9 )
    {
      cout << "  Enter the parameter ALPHA for this rule:\n";
      cin >> alpha_1d;
    }
    else
    {
      alpha_1d = 0.0;
    }

    if ( rule_1d == 9 )
    {
      cout << "  Enter the parameter BETA for this rule:\n";
      cin >> beta_1d;
    }
    else
    {
      beta_1d = 0.0;
    }

    for ( dim = dim_index; dim < dim_index + dim_inc; dim++ )
    {
      rule[dim] = rule_1d;
      alpha[dim] = alpha_1d;
      beta[dim] = beta_1d;
    }
    dim_index = dim_index + dim_inc;
  }
//
//  Get the filename.
//
  cout << "\n";
  cout << "  Enter an identifier to use for the filenames:\n";
  cin >> file_name;
//
//  Create the dataset.
//
  tol = sqrt ( webbur::r8_epsilon ( ) );
  sparse_grid_mixed_dataset_handle ( dim_num, level_max, rule, alpha, beta, 
    tol, file_name );
//
//  Terminate.
//
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  cout << "\n";
  cout << "SPARSE_GRID_MIXED_DATASET\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  webbur::timestamp ( );

  return 0;
}
//****************************************************************************80

char ch_cap ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= ch && ch <= 122 ) 
  {
    ch = ch - 32;
  }   

  return ch;
}
//***************************************************************************80

int rule_string_to_index ( string rule_string )

//***************************************************************************80
//
//  Purpose:
//
//    RULE_STRING_TO_INDEX converts a string identifying a rule to an index.
//
//  Discussion:
//
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
//     2, "F2",  Fejer Type 2, Open Fully Nested rule.
//     3, "GP",  Gauss Patterson, Open Fully Nested rule.
//     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
//     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
//     7, "LG",  Gauss Laguerre, Open Non Nested rule.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
//     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
//    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
//    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
//    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
//    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string RULE_STRING, a string.
//
//    Output, int RULE_STRING_TO_INDEX, the rule index.
//
{
  int rule_1d;

  if ( s_eqi ( rule_string, "CC" ) )
  {
    rule_1d = 1;
  }
  else if ( s_eqi ( rule_string, "F2" ) )
  {
    rule_1d = 2;
  }
  else if ( s_eqi ( rule_string, "GP" ) )
  {
    rule_1d = 3;
  }
  else if ( s_eqi ( rule_string, "GL" ) )
  {
    rule_1d = 4;
  }
  else if ( s_eqi ( rule_string, "GH" ) )
  {
    rule_1d = 5;
  }
  else if ( s_eqi ( rule_string, "GGH" ) )
  {
    rule_1d = 6;
  }
  else if ( s_eqi ( rule_string, "LG" ) )
  {
    rule_1d = 7;
  }
  else if ( s_eqi ( rule_string, "GLG" ) )
  {
    rule_1d = 8;
  }
  else if ( s_eqi ( rule_string, "GJ" ) )
  {
    rule_1d = 9;
  }
  else if ( s_eqi ( rule_string, "GW" ) )
  {
    rule_1d = 10;
  }
  else if ( s_eqi ( rule_string, "CCS" ) )
  {
    rule_1d = 11;
  }
  else if ( s_eqi ( rule_string, "F2S" ) )
  {
    rule_1d = 12;
  }
  else if ( s_eqi ( rule_string, "GPS" ) )
  {
    rule_1d = 13;
  }
  else
  {
    cout << "\n";
    cout << "RULE_STRING_TO_INDEX - Fatal error!\n";
    cout << "  Unexepected string.\n";
    exit ( 1 );
  }

  return rule_1d;
}
//****************************************************************************80

bool s_eqi ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
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
//  Parameters:
//
//    Input, string S1, S2, two strings.
//
//    Output, bool S_EQI, is true if the strings are equal. 
//
{
  int i;
  int nchar;
  int s1_length;
  int s2_length;

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  if ( s1_length < s2_length )
  {
    nchar = s1_length;
  }
  else
  {
    nchar = s2_length;
  }
//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ ) 
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) ) 
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < s1_length ) 
  {
    for ( i = nchar; i < s1_length; i++ ) 
    {
      if ( s1[i] != ' ' ) 
      {
        return false;
      }
    } 
  }
  else if ( nchar < s2_length ) 
  {
    for ( i = nchar; i < s2_length; i++ )
    {
      if ( s2[i] != ' ' ) 
      {
        return false;
      }
    } 
  }

  return true;
}
//***************************************************************************80

void sparse_grid_mixed_dataset_handle ( int dim_num, int level_max, int rule[], 
  double alpha[], double beta[], double tol, std::string file_name )

//***************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_MIXED_DATASET_HANDLE creates the dataset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer DIM_NUM, the spatial dimension.
//
//    Input, integer LEVEL_MAX, the level that defines the grid.
//
//    Input, int RULE[DIM_NUM], the rule in each dimension.
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
//     2, "F2",  Fejer Type 2, Open Fully Nested rule.
//     3, "GP",  Gauss Patterson, Open Fully Nested rule.
//     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
//     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
//     7, "LG",  Gauss Laguerre, Open Non Nested rule.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
//     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
//    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
//    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
//    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
//    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.
//
//    Input, double ALPHA[DIM_NUM], BETA[DIM_NUM], parameters used for
//    Generalized Gauss Hermite, Generalized Gauss Laguerre, 
//    and Gauss Jacobi rules.
//
//    Input, double TOL, a tolerance for point equality.
//
//    Input, string FILE_NAME, the main name of the 
//    output files.
//
{
  int point_num;
  int point_total_num;
  int *sparse_index;
  int *sparse_order;
  double *sparse_point;
  int *sparse_unique_index;
  double *sparse_weight;
//
//  Compute necessary data.
//
  point_total_num = webbur::sparse_grid_mixed_size_total ( dim_num, level_max, rule );

  point_num = webbur::sparse_grid_mixed_size ( dim_num, level_max, rule, alpha, 
    beta, tol );

  sparse_unique_index = new int[point_total_num];

  webbur::sparse_grid_mixed_unique_index ( dim_num, level_max, rule, alpha, beta,
    tol, point_num, point_total_num, sparse_unique_index );

  sparse_order = new int[dim_num*point_num];
  sparse_index = new int[dim_num*point_num];

  webbur::sparse_grid_mixed_index ( dim_num, level_max, rule, point_num, 
    point_total_num, sparse_unique_index, sparse_order, sparse_index );
//
//  Compute points and weights.
//
  sparse_point = new double [ dim_num * point_num ];

  webbur::sparse_grid_mixed_point ( dim_num, level_max, rule, alpha, beta, 
    point_num, sparse_order, sparse_index, sparse_point );

  sparse_weight = new double[point_num];

  webbur::sparse_grid_mixed_weight ( dim_num, level_max, rule, alpha, beta, 
    point_num, point_total_num, sparse_unique_index, sparse_weight );
//
//  Write points and weights to files.
//
  webbur::sparse_grid_mixed_write ( dim_num, rule, alpha, beta, point_num, 
    sparse_weight, sparse_point, file_name );

  delete [] sparse_index;
  delete [] sparse_order;
  delete [] sparse_point;
  delete [] sparse_unique_index;
  delete [] sparse_weight;

  return;
}
