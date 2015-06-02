# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# define UNBURNT 0
# define SMOLDERING 1
# define BURNING 2
# define BURNT 3

int main ( );
bool fire_spreads ( double prob_spread ) ;
int **forest_allocate ( int forest_size );
void forest_burns ( int forest_size, int **forest, double prob_spread );
void forest_delete ( int forest_size, int **forest );
void forest_initialize ( int forest_size, int **forest );
bool forest_is_burning ( int forest_size, int **forest );
void forest_print ( int forest_size, int **forest, int i_ignite, int j_ignite );
double get_percent_burned ( int forest_size, int **forest );
int seed_by_time ( int offset );
void timestamp ( void );
void tree_ignite ( int forest_size, int **forest, int i_ignite, 
  int j_ignite );

//****************************************************************************80

int main ( ) 

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FIRE_SERIAL.
//
//  Modified:
//
//    30 June 2013.
//
{
  int **forest;
  int forest_size = 20;
  int i;
  int i_ignite;
  int j_ignite;
  int offset;
  double percent_burned = 0.0;
  double prob_spread = 0.5;
  int seed;
  double u;

  timestamp ( );
  cout << "\n";
  cout << "FIRE_SERIAL\n";
  cout << "  C++ version\n";
  cout << "  A probabilistic simulation of a forest fire.\n";
  cout << "  The probability of tree-to-tree spread is " << prob_spread << "\n";
//
//  Initialize the random number generator.
//
  offset = 0;
  seed = seed_by_time ( offset );
  cout << "  The random number generator is seeded by " << seed << ".\n";
  srand ( seed );
//
//  Create the forest.
//
  forest = forest_allocate ( forest_size );
//
//  Initialize the values in the forest.
//
  forest_initialize ( forest_size, forest );
//
//  Choose a tree at random where the fire will start.
//
  u = ( double ) rand ( ) / ( double ) RAND_MAX;
  i_ignite = ( int ) ( ( double ) forest_size * u );
  u = ( double ) rand ( ) / ( double ) RAND_MAX;
  j_ignite = ( int ) ( ( double ) forest_size * u );
  tree_ignite ( forest_size, forest, i_ignite, j_ignite );
  cout << "\n";
  cout << "  Fire starts at tree[" << i_ignite << "][" << j_ignite << "].\n";
//
//  Let time run until nothing is burning any more.
//
  while ( forest_is_burning ( forest_size, forest ) )
  {
    forest_burns ( forest_size, forest, prob_spread );
  }
//
//  Display the final forest state.
//
  forest_print ( forest_size, forest, i_ignite, j_ignite );
//
//  Report the percentage of forest burned.
//
  percent_burned = get_percent_burned ( forest_size, forest );

  cout << "\n";
  cout << "  Percentage of forest burned = " << percent_burned << "\n";
//
//  Free memory.
//
  forest_delete ( forest_size, forest );
//
//  Terminate.
//
  cout << "\n";
  cout << "FIRE_SERIAL:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

bool fire_spreads ( double prob_spread ) 

//****************************************************************************80
//
//  Purpose:
//
//    FIRE_SPREADS determines whether the fire spreads.
//
//  Modified:
//
//    30 June 2013
//
//  Parameters:
//
//    Input, double PROB_SPREAD, the probability of spreading.
//
//    Output, bool FIRE_SPREADS, is TRUE if the fire spreads.
//
{
  double u;
  bool value;

  u = ( double ) rand ( ) / ( double ) RAND_MAX;

  if ( u < prob_spread )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

int **forest_allocate ( int forest_size )

//****************************************************************************80
//
//  Purpose:
//
//    FOREST_ALLOCATE allocates memory for a forest of the given size.
//
//  Modified:
//
//    30 June 2013
//
//  Parameters:
//
//    Input, int FOREST_SIZE, the linear dimension of the forest.
//
//    Output, int **FOREST_ALLOCATE, a pointer to the FOREST_SIZE x FOREST_SIZE
//    forest.
//
{
  int **forest;
  int i;
  int j;

  forest = new int *[forest_size];
  for ( i = 0; i < forest_size; i++ )
  {
    forest[i] = new int[forest_size];
  }

  return forest;
}
//****************************************************************************80

void forest_burns ( int forest_size, int **forest, double prob_spread )

//****************************************************************************80
//
//  Purpose:
//
//    FOREST_BURNS models a single time step of the burning forest.
//
//  Modified:
//
//    30 June 2013
//
//  Parameters:
//
//    Input, int FOREST_SIZE, the linear dimension of the forest.
//
//    Input, int **FOREST, a pointer to the FOREST_SIZE x FOREST_SIZE
//    forest.
//
//    Input, double PROB_SPREAD, the probability that the fire will spread
//    from a burning tree to an unburnt one.
//
{
  int i;
  int j;
//
//  Burning trees burn down;
//  Smoldering trees ignite;
//
  for ( i = 0; i < forest_size; i++ )
  {
    for ( j = 0; j < forest_size; j++ )
    {
      if ( forest[i][j] == BURNING )
      {
        forest[i][j] = BURNT;
      }
      else if ( forest[i][j] == SMOLDERING )
      {
        forest[i][j] = BURNING;
      }
    }
  }
//
//  Unburnt trees might catch fire.
//
  for ( i = 0; i < forest_size; i++ )
  {
    for ( j = 0; j < forest_size; j++ )
    {
      if ( forest[i][j] == BURNING)
      {
//
//  North.
//
        if ( i != 0 ) 
        {
          if ( fire_spreads ( prob_spread ) && forest[i-1][j] == UNBURNT )
          {
            forest[i-1][j] = SMOLDERING;
          }
        }
//
//  South.
//
        if ( i != forest_size - 1 ) 
        {
          if ( fire_spreads ( prob_spread ) && forest[i+1][j] == UNBURNT )
          {
            forest[i+1][j] = SMOLDERING;
          }
        }
//
//  West.
//
        if ( j != 0 )
        {
          if ( fire_spreads ( prob_spread ) && forest[i][j-1] == UNBURNT )
          {
            forest[i][j-1] = SMOLDERING;
          }
        }
//
//  East.
//
        if ( j != forest_size - 1 )
        {
          if ( fire_spreads ( prob_spread ) && forest[i][j+1] == UNBURNT )
          {
            forest[i][j+1] = SMOLDERING;
          }
        }

      }
    }
  }
  return;
}
//****************************************************************************80

void forest_delete ( int forest_size, int **forest ) 

//****************************************************************************80
//
//  Purpose:
//
//    FOREST_DELETE deletes the memory previously allocated for a forest.
//
//  Modified:
//
//    30 June 2013
//
//  Parameters:
//
//    Input, int FOREST_SIZE, the linear dimension of the forest.
//
//    Input, int **FOREST, a pointer to the FOREST_SIZE x FOREST_SIZE
//    forest.
//
{
  int i;

  for ( i = 0; i < forest_size; i++ )
  {
    delete [] forest[i];
  }
  delete [] forest;

  return;
}
//****************************************************************************80

void forest_initialize ( int forest_size, int ** forest ) 

//****************************************************************************80
//
//  Purpose:
//
//    FOREST_INITIALIZE initializes the forest values.
//
//  Modified:
//
//    30 June 2013
//
//  Parameters:
//
//    Input, int FOREST_SIZE, the linear dimension of the forest.
//
//    Input, int **FOREST, a pointer to the FOREST_SIZE x FOREST_SIZE
//    forest.  The entries in the forest have been initialized to "UNBURNT".
//
{
  int i;
  int j;

  for ( i = 0; i < forest_size; i++ )
  {
    for ( j = 0; j < forest_size; j++ )
    {
      forest[i][j] = UNBURNT;
    }
  }
  return;
}
//****************************************************************************80

bool forest_is_burning ( int forest_size, int **forest ) 

//****************************************************************************80
//
//  Purpose:
//
//    FOREST_IS_BURNING reports whether any trees in the forest are burning.
//
//  Modified:
//
//    30 June 2013
//
//  Parameters:
//
//    Input, int FOREST_SIZE, the linear dimension of the forest.
//
//    Input, int **FOREST, a pointer to the FOREST_SIZE x FOREST_SIZE
//    forest.
//
//    Output, bool FOREST_IS_BURNING, is TRUE if any tree in the forest
//    is in the SMOLDERING or BURNING state.
//
{
  int i;
  int j;
  bool value;

  value = false;

  for ( i = 0; i < forest_size; i++ ) 
  {
    for ( j = 0; j < forest_size; j++ )
    {
      if ( forest[i][j] == SMOLDERING || forest[i][j] == BURNING )
      {
        value = true;
        return value;
      }
    }
  }

  return value;
}
//****************************************************************************80

void forest_print ( int forest_size, int **forest, int i_ignite, int j_ignite )

//****************************************************************************80
//
//  Purpose:
//
//    FOREST_PRINT prints the state of the trees in the forest.
//
//  Modified:
//
//    30 June 2013
//
//  Parameters:
//
//    Input, int FOREST_SIZE, the linear dimension of the forest.
//
//    Input, int **FOREST, a pointer to the FOREST_SIZE x FOREST_SIZE
//    forest.
//
//    Input, int I_IGNITE, J_IGNITE, the location of the start of the fire.
//
{
  int i;
  int j;

  cout << "\n";
  cout << "  Map of fire damage.\n";
  cout << "  Fire started at '*'.\n";
  cout << "  Burned trees are indicated by '.'.\n";
  cout << "  Unburned trees are indicated by 'X'.\n";
  cout << "\n";

  for ( i = 0; i < forest_size; i++ )
  {
    cout << "  ";
    for ( j = 0; j < forest_size; j++ )
    {
      if ( i == i_ignite && j == j_ignite )
      {
        cout << "*";
      }
      else if ( forest[i][j] == BURNT )
      {
        cout << ".";
      } 
      else
      {
        cout << "X";
      }
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

double get_percent_burned ( int forest_size, int **forest ) 

//****************************************************************************80
//
//  Purpose:
//
//    GET_PERCENT_BURNED computes the percentage of the forest that burned.
//
//  Modified:
//
//    30 June 2013
//
//  Parameters:
//
//    Input, int FOREST_SIZE, the linear dimension of the forest.
//
//    Input, int **FOREST, a pointer to the FOREST_SIZE x FOREST_SIZE
//    forest.
//
//    Output, double GET_PERCENT_BURNED, the percentage of the forest
//    that burned.
//
{
  int i;
  int j;
  int sum;
  double value;

  sum = 0;
  for ( i = 0; i < forest_size; i++ ) 
  {
    for ( j = 0; j < forest_size; j++ )
    {
      if ( forest[i][j] == BURNT )
      {
        sum = sum + 1;
      }
    }
  }

  value = ( double ) ( sum ) / ( double ) ( forest_size * forest_size );

  return value;
}
//****************************************************************************80

int seed_by_time ( int offset ) 

//****************************************************************************80
//
//  Purpose:
//
//    SEED_BY_TIME seeds the random number generator using the time as a seed.
//
//  Modified:
//
//    30 June 2013;
//
//  Parameters:
//
//    Input, int OFFSET, an offset to add to the time.
//
//    Output, int SEED_BY_TIME, the seed used.
//
{
  int seed;
  time_t the_time;

  time ( &the_time );

  seed = ( int ) the_time + offset;

  return seed;
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
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void tree_ignite ( int forest_size, int **forest, int i_ignite, int j_ignite )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_IGNITE sets a given tree to the SMOLDERING state.
//
//  Modified:
//
//    30 June 2013
//
//  Parameters:
//
//    Input, int FOREST_SIZE, the linear dimension of the forest.
//
//    Input, int **FOREST, a pointer to the FOREST_SIZE x FOREST_SIZE
//    forest.
//
//    Input, int I_IGNITE, J_IGNITE, the coordinates of the tree which is
//    to be set to SMOLDERING.
//
{
  forest[i_ignite][j_ignite] = SMOLDERING;

  return;
}
